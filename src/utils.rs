use std::ops::Neg;

pub(crate) trait BucketIndex {
    type Out;
    fn get(window_index: usize, window_size: usize, bytes: &[u8]) -> Self::Out;
}

pub(crate) struct Booth;

pub(crate) struct Slice;

impl BucketIndex for Booth {
    type Out = i32;
    fn get(window_index: usize, window_size: usize, bytes: &[u8]) -> i32 {
        // Booth encoding:
        // * step by `window` size
        // * slice by size of `window + 1``
        // * each window overlap by 1 bit
        // * append a zero bit to the least significant end
        // Indexing rule for example window size 3 where we slice by 4 bits:
        // `[0, +1, +1, +2, +2, +3, +3, +4, -4, -3, -3 -2, -2, -1, -1, 0]``
        // So we can reduce the bucket size without preprocessing scalars
        // and remembering them as in classic signed digit encoding

        let skip_bits = (window_index * window_size).saturating_sub(1);
        let skip_bytes = skip_bits / 8;

        // fill into a u32
        let mut v: [u8; 4] = [0; 4];
        for (dst, src) in v.iter_mut().zip(bytes.iter().skip(skip_bytes)) {
            *dst = *src
        }
        let mut tmp = u32::from_le_bytes(v);

        // pad with one 0 if slicing the least significant window
        if window_index == 0 {
            tmp <<= 1;
        }

        // remove further bits
        tmp >>= skip_bits - (skip_bytes * 8);
        // apply the booth window
        tmp &= (1 << (window_size + 1)) - 1;

        let sign = tmp & (1 << window_size) == 0;

        // div ceil by 2
        tmp = (tmp + 1) >> 1;

        // find the booth action index
        if sign {
            tmp as i32
        } else {
            ((!(tmp - 1) & ((1 << window_size) - 1)) as i32).neg()
        }
    }
}

impl BucketIndex for Slice {
    type Out = usize;
    fn get(segment: usize, c: usize, bytes: &[u8]) -> usize {
        let skip_bits = segment * c;
        let skip_bytes = skip_bits / 8;

        if skip_bytes >= 32 {
            return 0;
        }

        let mut v = [0; 8];
        for (v, o) in v.iter_mut().zip(bytes.as_ref()[skip_bytes..].iter()) {
            *v = *o;
        }

        let mut tmp = u64::from_le_bytes(v);
        tmp >>= skip_bits - (skip_bytes * 8);
        tmp %= 1 << c;

        tmp as usize
    }
}
