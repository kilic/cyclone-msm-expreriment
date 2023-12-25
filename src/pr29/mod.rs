use ff::PrimeField;
use group::Group;
use halo2curves::bn256::{Fr, G1Affine, G1};

use self::round::Round;

#[macro_export]
macro_rules! div_ceil {
    ($a:expr, $b:expr) => {
        (($a - 1) / $b) + 1
    };
}

#[macro_export]
macro_rules! double_n {
    ($acc:expr, $n:expr) => {
        (0..$n).fold($acc, |acc, _| acc.double())
    };
}

#[macro_export]
macro_rules! range {
    ($index:expr, $n_items:expr) => {
        $index * $n_items..($index + 1) * $n_items
    };
}

#[macro_export]
macro_rules! index {
    ($digit:expr) => {
        ($digit & 0x7fffffff) as usize
    };
}

#[macro_export]
macro_rules! is_neg {
    ($digit:expr) => {
        sign_bit!($digit) != 0
    };
}

#[macro_export]
macro_rules! sign_bit {
    ($digit:expr) => {
        $digit & 0x80000000
    };
}

mod round;

pub struct MSM {
    signed_digits: Vec<u32>,
    sorted_positions: Vec<u32>,
    bucket_sizes: Vec<usize>,
    bucket_offsets: Vec<usize>,
    n_windows: usize,
    window: usize,
    n_buckets: usize,
    n_points: usize,
    round: Round,
}

impl MSM {
    pub fn allocate(n_points: usize, override_window: Option<usize>) -> Self {
        fn best_window(n: usize) -> usize {
            if n >= 262144 {
                15
            } else if n >= 65536 {
                12
            } else if n >= 16384 {
                11
            } else if n >= 8192 {
                10
            } else if n >= 1024 {
                9
            } else {
                7
            }
        }
        let window = match override_window {
            Some(window) => {
                let overriden = best_window(n_points);
                println!("override window from {} to {}", overriden, window);
                window
            }
            None => best_window(n_points),
        };
        let n_windows = div_ceil!(Fr::NUM_BITS as usize, window);
        let n_buckets = (1 << (window - 1)) + 1;
        let round = Round::new(n_buckets, n_points);
        MSM {
            signed_digits: vec![0u32; n_windows * n_points],
            sorted_positions: vec![0u32; n_windows * n_points],
            bucket_sizes: vec![0usize; n_windows * n_buckets],
            bucket_offsets: vec![0; n_buckets],
            n_windows,
            window,
            n_buckets,
            n_points,
            round,
        }
    }

    fn decompose(&mut self, scalars: &[Fr]) {
        pub(crate) fn get_bits(segment: usize, window: usize, bytes: &[u8]) -> u32 {
            let skip_bits = segment * window;
            let skip_bytes = skip_bits / 8;
            if skip_bytes >= 32 {
                return 0;
            }
            let mut v = [0; 4];
            for (v, o) in v.iter_mut().zip(bytes[skip_bytes..].iter()) {
                *v = *o;
            }
            let mut tmp = u32::from_le_bytes(v);
            tmp >>= skip_bits - (skip_bytes * 8);
            tmp %= 1 << window;
            tmp
        }
        let max = 1 << (self.window - 1);
        for (point_idx, scalar) in scalars.iter().enumerate() {
            let repr = scalar.to_repr();
            let mut borrow = 0u32;
            for window_idx in 0..self.n_windows {
                let windowed_digit = get_bits(window_idx, self.window, repr.as_ref()) + borrow;
                let signed_digit = if windowed_digit >= max {
                    borrow = 1;
                    ((1 << self.window) - windowed_digit) | 0x80000000
                } else {
                    borrow = 0;
                    windowed_digit
                };
                self.bucket_sizes[window_idx * self.n_buckets + index!(signed_digit)] += 1;
                self.signed_digits[window_idx * self.n_points + point_idx] = signed_digit;
            }
        }
        self.sort();
    }

    fn sort(&mut self) {
        // let t0 = start_timer!(|| "sort");
        for w_i in 0..self.n_windows {
            let sorted_positions = &mut self.sorted_positions[range!(w_i, self.n_points)];
            let bucket_sizes = &self.bucket_sizes[range!(w_i, self.n_buckets)];
            let signed_digits = &self.signed_digits[range!(w_i, self.n_points)];
            let mut offset = 0;
            for (i, size) in bucket_sizes.iter().enumerate() {
                self.bucket_offsets[i] = offset;
                offset += size;
            }
            for (sorted_idx, signed_digit) in signed_digits.iter().enumerate() {
                let bucket_idx = index!(signed_digit);
                sorted_positions[self.bucket_offsets[bucket_idx]] =
                    sign_bit!(signed_digit) | (sorted_idx as u32);
                self.bucket_offsets[bucket_idx] += 1;
            }
        }
        // end_timer!(t0);
    }

    pub fn evaluate_with(
        scalars: &[Fr],
        points: &[G1Affine],
        acc: &mut G1,
        override_window: Option<usize>,
    ) {
        let mut msm = Self::allocate(points.len(), override_window);
        // let t0 = start_timer!(|| "init");
        msm.decompose(scalars);
        // end_timer!(t0);
        for w_i in (0..msm.n_windows).rev() {
            if w_i != msm.n_windows - 1 {
                *acc = double_n!(*acc, msm.window);
            }
            // let t0 = start_timer!(|| "init");
            msm.round.init(
                points,
                &msm.sorted_positions[range!(w_i, msm.n_points)],
                &msm.bucket_sizes[range!(w_i, msm.n_buckets)],
            );
            // end_timer!(t0);
            // let t0 = start_timer!(|| "buckets");

            let buckets = msm.round.evaluate();
            // end_timer!(t0);
            // println!("buckets: {:?}", buckets.len());
            // let t0 = start_timer!(|| "sum");
            let mut running_sum = G1::identity();
            for bucket in buckets.into_iter().skip(1).rev() {
                running_sum += bucket;
                *acc += &running_sum;
            }
            // end_timer!(t0);
        }
    }
}
