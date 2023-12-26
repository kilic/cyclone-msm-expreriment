use core::panic;
use ff::{Field, PrimeField};
use group::{cofactor::CofactorCurveAffine, Group};
use halo2curves::{
    bn256::{Fq, Fr, G1Affine, G1},
    CurveAffine,
};
use rayon::{
    current_num_threads,
    iter::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator},
    scope,
};
use std::ops::Neg;

use crate::utils::{Booth, BucketIndex};

fn batch_add(
    size: usize,
    buckets: &mut [Bucket<G1Affine>],
    points: &[SchedulePoint],
    bases: &[G1Affine],
) {
    let mut t = vec![Fq::ZERO; size];
    let mut z = vec![Fq::ZERO; size];
    let mut acc = Fq::ONE;

    for (
        (
            SchedulePoint {
                base_idx,
                buck_idx,
                sign,
            },
            t,
        ),
        z,
    ) in points.iter().zip(t.iter_mut()).zip(z.iter_mut())
    {
        *z = buckets[*buck_idx].x() - bases[*base_idx].x;
        if *sign {
            *t = acc * (buckets[*buck_idx].y() - bases[*base_idx].y);
        } else {
            *t = acc * (buckets[*buck_idx].y() + bases[*base_idx].y);
        }
        acc *= *z;
    }

    acc = acc.invert().expect(":(");

    for (
        (
            SchedulePoint {
                base_idx,
                buck_idx,
                sign,
            },
            t,
        ),
        z,
    ) in points.iter().zip(t.iter()).zip(z.iter()).rev()
    {
        let lambda = acc * t;
        acc *= z;

        let x = lambda.square() - (buckets[*buck_idx].x() + bases[*base_idx].x);
        if *sign {
            buckets[*buck_idx].set_y(&((lambda * (bases[*base_idx].x - x)) - bases[*base_idx].y));
        } else {
            buckets[*buck_idx].set_y(&((lambda * (bases[*base_idx].x - x)) + bases[*base_idx].y));
        }
        buckets[*buck_idx].set_x(&x);
    }
}

#[derive(Debug, Clone)]
enum Bucket<T> {
    None,
    Point(T),
}

impl Bucket<G1> {
    fn add_assign(&mut self, point: &G1Affine, sign: bool) {
        *self = match *self {
            Bucket::None => Bucket::Point({
                if sign {
                    point.to_curve()
                } else {
                    point.to_curve().neg()
                }
            }),
            Bucket::Point(a) => {
                if sign {
                    Self::Point(a + point)
                } else {
                    Self::Point(a - point)
                }
            }
        }
    }

    fn add(&self, other: &Bucket<G1Affine>) -> G1 {
        match (self, other) {
            (Bucket::Point(this), Bucket::Point(other)) => this + other,
            (Bucket::Point(this), Bucket::None) => *this,
            (Bucket::None, Bucket::Point(other)) => other.to_curve(),
            (Bucket::None, Bucket::None) => G1::identity(),
        }
    }
}

impl Bucket<G1Affine> {
    fn assign(&mut self, point: &G1Affine, sign: bool) -> bool {
        match *self {
            Bucket::None => {
                *self = Bucket::Point(if sign { *point } else { point.neg() });
                true
            }
            Bucket::Point(_) => false,
        }
    }

    fn x(&self) -> Fq {
        match *self {
            Bucket::None => panic!("::x None"),
            Bucket::Point(a) => a.x,
        }
    }

    fn y(&self) -> Fq {
        match *self {
            Bucket::None => panic!("::y None"),
            Bucket::Point(a) => a.y,
        }
    }

    fn set_x(&mut self, x: &Fq) {
        match *self {
            Bucket::None => panic!("::set_x None"),
            Bucket::Point(ref mut a) => a.x = *x,
        }
    }

    fn _is_on_curve(&self) {
        match *self {
            Bucket::None => panic!("::_is_on_curve None"),
            Bucket::Point(a) => assert!(bool::from(a.is_on_curve())),
        }
    }

    fn set_y(&mut self, y: &Fq) {
        match *self {
            Bucket::None => panic!("::set_y None"),
            Bucket::Point(ref mut a) => a.y = *y,
        }
    }
}

struct Schedule {
    buckets: Vec<Bucket<G1Affine>>,
    set: Vec<SchedulePoint>,
    ptr: usize,
}

#[derive(Debug, Clone, Default)]
struct SchedulePoint {
    base_idx: usize,
    buck_idx: usize,
    sign: bool,
}

impl SchedulePoint {
    fn new(base_idx: usize, buck_idx: usize, sign: bool) -> Self {
        Self {
            base_idx,
            buck_idx,
            sign,
        }
    }
}

impl Schedule {
    fn new(batch_size: usize, c: usize) -> Self {
        Self {
            buckets: vec![Bucket::None; 1 << (c - 1)],
            set: vec![SchedulePoint::default(); batch_size],
            ptr: 0,
        }
    }

    fn contains(&self, buck_idx: usize) -> bool {
        // Note:
        // * BTreeMap is slower
        // * HashMap  is not a double ended iterator
        self.set
            .iter()
            .position(|sch| sch.buck_idx == buck_idx)
            .is_some()
    }

    fn execute(&mut self, bases: &[G1Affine]) {
        if self.ptr != 0 {
            batch_add(self.ptr, &mut self.buckets, &self.set, bases);
            self.ptr = 0;
            self.set
                .iter_mut()
                .for_each(|sch| *sch = SchedulePoint::default());
        }
    }

    fn add(&mut self, bases: &[G1Affine], base_idx: usize, buck_idx: usize, sign: bool) {
        if !self.buckets[buck_idx].assign(&bases[base_idx], sign) {
            self.set[self.ptr] = SchedulePoint::new(base_idx, buck_idx, sign);
            self.ptr += 1;
        }

        if self.ptr == self.set.len() {
            self.execute(bases);
        }
    }
}

pub fn msm_serial(coeffs: &[Fr], bases: &[G1Affine], c: usize, batch_size: usize) -> G1 {
    let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();

    let segments = (256 / c) + 1;
    let mut acc = G1::identity();

    for current_segment in (0..segments).rev() {
        for _ in 0..c {
            acc = acc.double();
        }

        let mut j_bucks = vec![Bucket::None; 1 << (c - 1)];
        let mut sched = Schedule::new(batch_size, c);

        for (base_idx, coeff) in coeffs.iter().enumerate() {
            let buck_idx = Booth::get(current_segment, c, coeff.as_ref());

            if buck_idx != 0 {
                let sign = buck_idx.is_positive();
                let buck_idx = buck_idx.unsigned_abs() as usize - 1;

                if sched.contains(buck_idx) {
                    j_bucks[buck_idx].add_assign(&bases[base_idx], sign);
                } else {
                    sched.add(bases, base_idx, buck_idx, sign);
                }
            }
        }

        sched.execute(bases);

        let mut running_sum = G1::identity();
        for (j_buck, a_buck) in j_bucks.iter().zip(sched.buckets.iter()).rev() {
            running_sum += j_buck.add(a_buck);
            acc += &running_sum;
        }
    }
    acc
}

pub fn msm_par_split(coeffs: &[Fr], bases: &[G1Affine], c: usize, batch_size: usize) -> G1 {
    assert_eq!(coeffs.len(), bases.len());

    let num_threads = current_num_threads();

    let chunk = coeffs.len() / num_threads;
    let num_chunks = coeffs.chunks(chunk).len();
    let mut results = vec![G1::identity(); num_chunks];
    scope(|scope| {
        let chunk = coeffs.len() / num_threads;

        for ((coeffs, bases), acc) in coeffs
            .chunks(chunk)
            .zip(bases.chunks(chunk))
            .zip(results.iter_mut())
        {
            scope.spawn(move |_| {
                *acc = msm_serial(coeffs, bases, c, batch_size);
            });
        }
    });
    results.iter().sum()
}

pub fn msm_par_window(coeffs: &[Fr], bases: &[G1Affine], c: usize, batch_size: usize) -> G1 {
    let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();

    let segments = (256 / c) + 1;
    let mut acc = vec![G1::identity(); segments];
    acc.par_iter_mut().enumerate().rev().for_each(|(w, acc)| {
        let mut j_bucks = vec![Bucket::None; 1 << (c - 1)];
        let mut sched = Schedule::new(batch_size, c);

        for (base_idx, coeff) in coeffs.iter().enumerate() {
            let buck_idx = Booth::get(w, c, coeff.as_ref());

            if buck_idx != 0 {
                let sign = buck_idx.is_positive();
                let buck_idx = buck_idx.unsigned_abs() as usize - 1;

                if sched.contains(buck_idx) {
                    j_bucks[buck_idx].add_assign(&bases[base_idx], sign);
                } else {
                    sched.add(bases, base_idx, buck_idx, sign);
                }
            }
        }

        sched.execute(bases);

        let mut running_sum = G1::identity();
        for (j_buck, a_buck) in j_bucks.iter().zip(sched.buckets.iter()).rev() {
            running_sum += j_buck.add(a_buck);
            *acc += running_sum;
        }
        for _ in 0..c * w {
            *acc = acc.double();
        }
    });
    acc.into_iter().sum::<G1>()
}
