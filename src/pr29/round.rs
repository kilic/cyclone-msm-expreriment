use ff::Field;
use group::cofactor::CofactorCurveAffine;
use halo2curves::bn256::{Fq, G1Affine};

macro_rules! log_ceil {
    ($a:expr) => {
        $a.next_power_of_two().trailing_zeros()
    };
}

pub(crate) struct Round {
    t: Vec<G1Affine>,
    odd_points: Vec<G1Affine>,
    bucket_sizes: Vec<usize>,
    n_buckets: usize,
    n_points: usize,
    out_off: usize,
    in_off: usize,
}

macro_rules! first_phase {
    ($out:expr, $p0:expr, $p1:expr, $acc:expr) => {
        $out.x = $p0.x + $p1.x;
        $p1.x = $p1.x - $p0.x;
        $p1.y = ($p1.y - $p0.y) * $acc;
        $acc = $acc * $p1.x;
    };
}
macro_rules! second_phase {
    ($out:expr, $p0:expr, $p1:expr, $acc:expr) => {
        $p1.y = $p1.y * $acc;
        $acc = $acc * $p1.x;
        $out.x = $p1.y.square() - $out.x;
        $out.y = ($p1.y * ($p0.x - $out.x)) - $p0.y;
    };
}

impl Round {
    pub(crate) fn new(n_buckets: usize, n_points: usize) -> Self {
        let odd_points = vec![G1Affine::identity(); n_buckets];
        let bucket_sizes = vec![0; n_buckets];
        //  TODO: requires less than allocated
        let t = vec![G1Affine::identity(); n_points * 2];
        Self {
            t,
            odd_points,
            bucket_sizes,
            n_buckets,
            n_points,
            out_off: 0,
            in_off: 0,
        }
    }
    pub(crate) fn init(&mut self, bases: &[G1Affine], positions: &[u32], bucket_sizes: &[usize]) {
        {
            assert_eq!(self.n_points, positions.len());
            assert_eq!(self.n_points, bases.len());
            assert_eq!(self.bucket_sizes.len(), self.n_buckets);
        }
        macro_rules! get_base {
            ($positions:expr, $off:expr) => {
                if is_neg!($positions[$off]) {
                    -bases[index!($positions[$off])]
                } else {
                    bases[index!($positions[$off])]
                }
            };
        }
        let n_additions = bucket_sizes
            .iter()
            .map(|bucket_sizes| bucket_sizes / 2)
            .collect::<Vec<usize>>();
        let mut out_off = n_additions.iter().sum::<usize>();
        let mut tmp_off = 0;
        let mut position_off = 0;
        let mut acc = Fq::ONE;

        for (bucket_index, bucket_size) in bucket_sizes.iter().enumerate() {
            let positions = &positions[position_off..position_off + bucket_size];
            position_off += bucket_size;
            let bucket_size_pre = positions.len();
            if bucket_size_pre == 0 {
                self.odd_points[bucket_index] = G1Affine::identity();
                self.bucket_sizes[bucket_index] = 0;
            } else {
                let mut in_off = 0;
                let bucket_size_post = (bucket_size_pre + 1) / 2;
                let n_additions = bucket_size_pre / 2;
                // t_inv.resize(10000, Fq::ZERO);
                // process even number of additions
                for _ in 0..n_additions & (usize::MAX - 1) {
                    // second operand must be mutable
                    self.t[tmp_off] = get_base!(positions, in_off);
                    let lhs = get_base!(positions, in_off + 1);
                    first_phase!(self.t[out_off], lhs, self.t[tmp_off], acc);
                    tmp_off += 1;
                    out_off += 1;
                    in_off += 2;
                }
                // process the latest elements if there are odd number of additions
                match (bucket_size_pre & 1 == 1, bucket_size_post & 1 == 1) {
                    // 1 base point left
                    // move to odd-point cache
                    (true, true) => {
                        assert_eq!(positions.len() - 1, in_off);
                        self.odd_points[bucket_index] = get_base!(positions, in_off);
                    }
                    // 2 base point left
                    // move addition result to odd-point cache
                    (false, true) => {
                        self.t[tmp_off] = get_base!(positions, in_off);
                        let lhs = get_base!(positions, in_off + 1);
                        first_phase!(self.odd_points[bucket_index], lhs, self.t[tmp_off], acc);
                        tmp_off += 1;
                    }
                    // 3 base point left
                    // move addition of first two to intermediate and last to odd-point cache
                    (true, false) => {
                        self.t[tmp_off] = get_base!(positions, in_off);
                        let lhs = get_base!(positions, in_off + 1);
                        first_phase!(self.t[out_off], lhs, self.t[tmp_off], acc);
                        self.t[out_off + 1] = get_base!(positions, in_off + 2);
                        tmp_off += 1;
                        out_off += 2;
                    }
                    _ => { /* 0 base point left */ }
                }
                self.bucket_sizes[bucket_index] = bucket_size_post;
            }
        }
        self.in_off = tmp_off;
        self.out_off = out_off;
        tmp_off -= 1;
        out_off -= 1;

        acc = acc.invert().unwrap();

        for (bucket_index, bucket_size_pre) in bucket_sizes.iter().enumerate().rev() {
            let positions = &positions[position_off - bucket_size_pre..position_off];
            position_off -= bucket_size_pre;
            if *bucket_size_pre == 0 {
                // already updated in first phase
            } else {
                if positions.len() != 0 {
                    let bucket_size_post = (bucket_size_pre + 1) / 2;
                    let n_additions = bucket_size_pre / 2;
                    let mut in_off = positions.len() - 1;
                    // process the latest elements if there are odd number of additions
                    match (bucket_size_pre & 1 == 1, bucket_size_post & 1 == 1) {
                        // 1 base point left
                        // move to odd-point cache
                        (true, true) => {
                            in_off -= 1;
                        }
                        // 2 base point left
                        // move addition result to odd-point cache
                        (false, true) => {
                            let lhs = get_base!(positions, in_off);
                            second_phase!(self.odd_points[bucket_index], lhs, self.t[tmp_off], acc);
                            tmp_off -= 1;
                            in_off -= 2;
                        }
                        // 3 base point left
                        // move addition of first two to intermediate and last to odd-point cache
                        (true, false) => {
                            in_off -= 1;
                            out_off -= 1;
                            let lhs = get_base!(positions, in_off);
                            second_phase!(self.t[out_off], lhs, self.t[tmp_off], acc);
                            tmp_off -= 1;
                            in_off -= 2;
                            out_off -= 1;
                        }
                        _ => { /* 0 base point left */ }
                    }
                    // process even number of additions
                    for _ in (0..n_additions & (usize::MAX - 1)).rev() {
                        let lhs = get_base!(positions, in_off);
                        second_phase!(self.t[out_off], lhs, self.t[tmp_off], acc);
                        tmp_off -= 1;
                        out_off -= 1;
                        in_off -= 2;
                    }
                }
            }
        }
    }
    fn batch_add(&mut self) {
        let (mut out_off, mut in_off) = (self.out_off, self.in_off);
        let mut acc = Fq::ONE;

        for bucket_index in 0..self.n_buckets {
            let bucket_size_pre = self.bucket_sizes[bucket_index];
            let n_additions = bucket_size_pre / 2;
            let bucket_size_post = (bucket_size_pre + 1) / 2;
            for _ in 0..n_additions & (usize::MAX - 1) {
                first_phase!(self.t[out_off], self.t[in_off], self.t[in_off + 1], acc);
                (out_off, in_off) = (out_off + 1, in_off + 2);
            }
            match (bucket_size_pre & 1 == 1, bucket_size_post & 1 == 1) {
                (true, false) => {
                    first_phase!(self.t[out_off], self.t[in_off], self.t[in_off + 1], acc);
                    self.t[out_off + 1] = self.odd_points[bucket_index];
                    out_off += 2;
                    in_off += 2;
                }
                (false, true) => {
                    first_phase!(
                        self.odd_points[bucket_index],
                        self.t[in_off],
                        self.t[in_off + 1],
                        acc
                    );
                    in_off += 2;
                }
                _ => { /* clean sheets */ }
            }
        }
        self.out_off = out_off;
        self.in_off = in_off;
        out_off -= 1;
        in_off -= 2;
        acc = acc.invert().unwrap();
        // process second phase
        for bucket_index in (0..self.n_buckets).rev() {
            let bucket_size_pre = self.bucket_sizes[bucket_index];
            let n_additions = bucket_size_pre / 2;
            let bucket_size_post = (bucket_size_pre + 1) / 2;

            match (bucket_size_pre & 1 == 1, bucket_size_post & 1 == 1) {
                (true, false) => {
                    out_off -= 1;
                    second_phase!(self.t[out_off], self.t[in_off], self.t[in_off + 1], acc);
                    out_off -= 1;
                    in_off -= 2;
                }
                (false, true) => {
                    second_phase!(
                        self.odd_points[bucket_index],
                        self.t[in_off],
                        self.t[in_off + 1],
                        acc
                    );
                    in_off -= 2;
                }
                _ => { /* clean sheets */ }
            }

            for _ in 0..n_additions & (usize::MAX - 1) {
                second_phase!(self.t[out_off], self.t[in_off], self.t[in_off + 1], acc);
                out_off -= 1;
                in_off -= 2;
            }
            self.bucket_sizes[bucket_index] = bucket_size_post;
        }
    }

    fn max_tree_height(&self) -> usize {
        *self.tree_heights().iter().max().unwrap()
    }
    fn tree_heights(&self) -> Vec<usize> {
        self.bucket_sizes
            .iter()
            .map(|bucket_size| log_ceil!(bucket_size) as usize)
            .collect()
    }
    pub fn evaluate(&mut self) -> &[G1Affine] {
        for _ in 0..self.max_tree_height() {
            self.batch_add();
        }
        &self.odd_points
    }
}
