use std::ops::Neg;

use ark_std::{end_timer, start_timer};
use ff::{Field, PrimeField};
use group::{cofactor::CofactorCurveAffine, Curve, Group};
use halo2curves::bn256::{Fq, Fr, G1Affine, G1};

use crate::utils::{get_booth_index, get_bucket_index};

#[derive(Debug, Clone, Default)]
struct Sched {
    base_idx: usize,
    buck_idx: usize,
    sign: bool,
}

impl Sched {
    fn contains(set: &[Self], buck_idx: usize) -> bool {
        set.iter()
            .position(|sched| sched.buck_idx == buck_idx)
            .is_some()
    }

    fn clear(set: &mut [Self]) {
        set.iter_mut().for_each(|sched| *sched = Sched::default());
    }

    fn batch_add(sched: &mut [Self], size: usize, buckets: &mut [G1Affine], bases: &[G1Affine]) {
        batch_add2(size, buckets, &sched, bases);
        Self::clear(sched);
    }

    fn new(base_idx: usize, buck_idx: usize, sign: bool) -> Self {
        Self {
            base_idx,
            buck_idx,
            sign,
        }
    }
}

fn batch_add(
    size: usize,
    buckets: &mut [G1Affine],
    sched_bases: &[(usize, usize)],
    bases: &[G1Affine],
) {
    let mut t = vec![Fq::ZERO; size];
    // TODO/Note: with z we removed one subtraction but is it worth allocating `size` mem?
    let mut z = vec![Fq::ZERO; size];
    let mut acc = Fq::ONE;

    for (((buck_idx, base_idx), t), z) in sched_bases.iter().zip(t.iter_mut()).zip(z.iter_mut()) {
        *z = buckets[*buck_idx].x - bases[*base_idx].x;
        *t = acc * (buckets[*buck_idx].y - bases[*base_idx].y);
        acc *= *z;
    }

    acc = acc.invert().unwrap();

    for (((buck_idx, base_idx), t), z) in sched_bases.iter().zip(t.iter()).zip(z.iter()).rev() {
        let lambda = acc * t;
        // let base = bases[*base_idx];
        acc *= z;

        let x = lambda.square() - (buckets[*buck_idx].x + bases[*base_idx].x);
        buckets[*buck_idx].y = (lambda * (bases[*base_idx].x - x)) - bases[*base_idx].y;
        buckets[*buck_idx].x = x;
    }
}

fn batch_add2(size: usize, buckets: &mut [G1Affine], sched_bases: &[Sched], bases: &[G1Affine]) {
    let mut t = vec![Fq::ZERO; size];
    // TODO/Note: with z we dedup one subtraction which is `p0.x - p1.x` but is it worth allocating `size` mem?
    let mut z = vec![Fq::ZERO; size];
    let mut acc = Fq::ONE;

    for (
        (
            Sched {
                base_idx,
                buck_idx,
                sign,
            },
            t,
        ),
        z,
    ) in sched_bases.iter().zip(t.iter_mut()).zip(z.iter_mut())
    {
        *z = buckets[*buck_idx].x - bases[*base_idx].x;
        if *sign {
            *t = acc * (buckets[*buck_idx].y - bases[*base_idx].y);
        } else {
            *t = acc * (buckets[*buck_idx].y + bases[*base_idx].y);
        }
        acc *= *z;
    }

    acc = acc.invert().unwrap();

    for (
        (
            Sched {
                base_idx,
                buck_idx,
                sign,
            },
            t,
        ),
        z,
    ) in sched_bases.iter().zip(t.iter()).zip(z.iter()).rev()
    {
        let lambda = acc * t;
        // let base = bases[*base_idx];
        acc *= z;

        let x = lambda.square() - (buckets[*buck_idx].x + bases[*base_idx].x);
        if *sign {
            buckets[*buck_idx].y = (lambda * (bases[*base_idx].x - x)) - bases[*base_idx].y;
        } else {
            buckets[*buck_idx].y = (lambda * (bases[*base_idx].x - x)) + bases[*base_idx].y;
        }
        buckets[*buck_idx].x = x;
    }
}

fn batch_add3(
    size: usize,
    buckets: &mut [Bucket<G1Affine>],
    sched_bases: &[Sched],
    bases: &[G1Affine],
) {
    let mut t = vec![Fq::ZERO; size];
    // TODO/Note: with z we dedup one subtraction which is `p0.x - p1.x` but is it worth allocating `size` mem?
    let mut z = vec![Fq::ZERO; size];
    let mut acc = Fq::ONE;

    for (
        (
            Sched {
                base_idx,
                buck_idx,
                sign,
            },
            t,
        ),
        z,
    ) in sched_bases.iter().zip(t.iter_mut()).zip(z.iter_mut())
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
            Sched {
                base_idx,
                buck_idx,
                sign,
            },
            t,
        ),
        z,
    ) in sched_bases.iter().zip(t.iter()).zip(z.iter()).rev()
    {
        let lambda = acc * t;
        // let base = bases[*base_idx];
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

pub fn multiexp_serial(coeffs: &[Fr], bases: &[G1Affine], c: usize, batch_size: usize) -> G1 {
    let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();

    let segments = (256 / c) + 1;
    let mut acc = G1::identity();

    let mut sched_bucks = vec![(0, 0); batch_size];
    let mut ptr = 0;

    for current_segment in (0..segments).rev() {
        // let mut stat = Stat::default();
        for _ in 0..c {
            acc = acc.double();
        }

        let mut a_bucks = vec![G1Affine::identity(); (1 << c) - 1];
        let mut j_bucks = vec![G1::identity(); (1 << c) - 1];

        for (base_idx, coeff) in coeffs.iter().enumerate() {
            let buck_idx = get_bucket_index(current_segment, c, coeff.as_ref());

            if buck_idx != 0 {
                let buck_idx = buck_idx - 1;

                if sched_bucks
                    .iter()
                    .position(|(sched_idx, _)| *sched_idx == buck_idx)
                    .is_some()
                {
                    // stat.n_jac += 1;
                    j_bucks[buck_idx] += bases[base_idx];
                } else {
                    // stat.n_aff += 1;

                    if bool::from(a_bucks[buck_idx].is_identity()) {
                        a_bucks[buck_idx] = bases[base_idx];
                    } else {
                        sched_bucks[ptr] = (buck_idx, base_idx);
                        ptr += 1;
                    }
                }
            }

            if ptr == batch_size {
                let mut t = vec![Fq::ZERO; batch_size];
                let mut acc = Fq::ONE;

                for ((buck_idx, base_idx), t) in sched_bucks.iter().zip(t.iter_mut()) {
                    *t = acc;
                    acc *= a_bucks[*buck_idx].x - bases[*base_idx].x;
                }

                acc = acc.invert().unwrap();

                for ((buck_idx, base_idx), t) in sched_bucks.iter().zip(t.iter()).rev() {
                    let inv = acc * t;
                    let base = bases[*base_idx];
                    acc = acc * (a_bucks[*buck_idx].x - base.x);

                    let lambda = (a_bucks[*buck_idx].y - base.y) * inv;
                    let x = lambda.square() - (a_bucks[*buck_idx].x + base.x);
                    a_bucks[*buck_idx].y = (lambda * (base.x - x)) - bases[*base_idx].y;
                    a_bucks[*buck_idx].x = x;
                }

                ptr = 0;
            }
        }

        {
            let mut t = vec![Fq::ZERO; ptr];
            let mut acc = Fq::ONE;

            for ((buck_idx, base_idx), t) in sched_bucks.iter().zip(t.iter_mut()).take(ptr) {
                *t = acc;
                acc *= a_bucks[*buck_idx].x - bases[*base_idx].x;
            }

            acc = acc.invert().unwrap();

            for ((buck_idx, base_idx), t) in sched_bucks.iter().zip(t.iter()).take(ptr).rev() {
                let inv = acc * t;
                let base = bases[*base_idx];
                acc = acc * (a_bucks[*buck_idx].x - base.x);

                let lambda = (a_bucks[*buck_idx].y - base.y) * inv;
                let x = lambda.square() - (a_bucks[*buck_idx].x + base.x);
                a_bucks[*buck_idx].y = (lambda * (base.x - x)) - bases[*base_idx].y;
                a_bucks[*buck_idx].x = x;
            }
            ptr = 0;
        }

        let mut running_sum = G1::identity();
        for (j_buck, a_buck) in j_bucks.into_iter().zip(a_bucks.into_iter()).rev() {
            running_sum = a_buck + j_buck + running_sum;
            acc += &running_sum;
        }
    }

    acc
}

pub fn multiexp_serial2(coeffs: &[Fr], bases: &[G1Affine], c: usize, batch_size: usize) -> G1 {
    let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();

    let segments = (256 / c) + 1;
    let mut acc = G1::identity();

    let mut sched_bases = vec![(0, 0); batch_size];
    let mut ptr = 0;

    for current_segment in (0..segments).rev() {
        for _ in 0..c {
            acc = acc.double();
        }

        let mut a_bucks = vec![G1Affine::identity(); (1 << c) - 1];
        let mut j_bucks = vec![G1::identity(); (1 << c) - 1];

        for (base_idx, coeff) in coeffs.iter().enumerate() {
            let buck_idx = get_bucket_index(current_segment, c, coeff.as_ref());

            if buck_idx != 0 {
                let buck_idx = buck_idx - 1;

                if sched_bases
                    .iter()
                    .position(|(sched_idx, _)| *sched_idx == buck_idx)
                    .is_some()
                {
                    j_bucks[buck_idx] += bases[base_idx];
                } else {
                    if bool::from(a_bucks[buck_idx].is_identity()) {
                        a_bucks[buck_idx] = bases[base_idx];
                    } else {
                        sched_bases[ptr] = (buck_idx, base_idx);
                        ptr += 1;
                    }
                }
            }

            if ptr == batch_size {
                // let t0 = start_timer!(|| "batch_add");
                batch_add(batch_size, &mut a_bucks, &sched_bases, bases);
                // end_timer!(t0);
                ptr = 0;
            }
        }

        batch_add(ptr, &mut a_bucks, &sched_bases, bases);
        ptr = 0;

        let mut running_sum = G1::identity();
        for (j_buck, a_buck) in j_bucks.into_iter().zip(a_bucks.into_iter()).rev() {
            running_sum += a_buck + j_buck;
            acc += &running_sum;
        }
    }

    acc
}

pub fn multiexp_serial3(coeffs: &[Fr], bases: &[G1Affine], c: usize, batch_size: usize) -> G1 {
    let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();

    let segments = (256 / c) + 1;
    let mut acc = G1::identity();

    let mut sched_bases = vec![Sched::default(); batch_size];
    let mut ptr = 0;

    // let mut stat = Stat::default();

    for current_segment in (0..segments).rev() {
        for _ in 0..c {
            acc = acc.double();
        }

        let mut a_bucks = vec![G1Affine::identity(); 1 << (c - 1)];
        let mut j_bucks = vec![G1::identity(); 1 << (c - 1)];

        for (base_idx, coeff) in coeffs.iter().enumerate() {
            let buck_idx = get_booth_index(current_segment, c, coeff.as_ref());

            if buck_idx != 0 {
                let sign = buck_idx.is_positive();
                let buck_idx = buck_idx.unsigned_abs() as usize - 1;

                if Sched::contains(&sched_bases, buck_idx as usize) {
                    // stat.add(current_segment, 0, 0, 1);
                    if sign {
                        j_bucks[buck_idx] += bases[base_idx];
                    } else {
                        j_bucks[buck_idx] -= bases[base_idx];
                    }
                } else {
                    // stat.add(current_segment, 0, 1, 0);
                    if bool::from(a_bucks[buck_idx].is_identity()) {
                        if sign {
                            a_bucks[buck_idx] = bases[base_idx];
                        } else {
                            a_bucks[buck_idx] = bases[base_idx].neg();
                        }
                    } else {
                        sched_bases[ptr] = Sched::new(base_idx, buck_idx, sign);
                        ptr += 1;
                    }
                }
            }

            if ptr == batch_size {
                // let t0 = start_timer!(|| "batch_add2");
                batch_add2(batch_size, &mut a_bucks, &sched_bases, bases);
                // end_timer!(t0);
                ptr = 0;
            }
        }

        {
            batch_add2(ptr, &mut a_bucks, &sched_bases, bases);
            ptr = 0;
        }

        let mut running_sum = G1::identity();
        for (j_buck, a_buck) in j_bucks.into_iter().zip(a_bucks.into_iter()).rev() {
            running_sum += a_buck + j_buck;
            acc += &running_sum;
        }
    }

    // stat.debug();
    acc
}

pub fn multiexp_serial4(coeffs: &[Fr], bases: &[G1Affine], c: usize, batch_size: usize) -> G1 {
    let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();

    let segments = (256 / c) + 1;
    let mut acc = G1::identity();

    let mut sched_bases = vec![Sched::default(); batch_size];
    let mut ptr = 0;

    // let mut stat = Stat::default();

    for current_segment in (0..segments).rev() {
        for _ in 0..c {
            acc = acc.double();
        }

        let mut a_bucks = vec![G1Affine::identity(); 1 << (c - 1)];
        let mut j_bucks = vec![G1::identity(); 1 << (c - 1)];

        for (base_idx, coeff) in coeffs.iter().enumerate() {
            let buck_idx = get_booth_index(current_segment, c, coeff.as_ref());

            if buck_idx != 0 {
                let sign = buck_idx.is_positive();
                let buck_idx = buck_idx.unsigned_abs() as usize - 1;

                if Sched::contains(&sched_bases, buck_idx as usize) {
                    // stat.add(current_segment, 0, 0, 1);
                    if sign {
                        j_bucks[buck_idx] += bases[base_idx];
                    } else {
                        j_bucks[buck_idx] += bases[base_idx].neg();
                    }
                } else {
                    // stat.add(current_segment, 0, 1, 0);
                    if bool::from(a_bucks[buck_idx].is_identity()) {
                        if sign {
                            a_bucks[buck_idx] = bases[base_idx];
                        } else {
                            a_bucks[buck_idx] = bases[base_idx].neg();
                        }
                    } else {
                        sched_bases[ptr] = Sched::new(base_idx, buck_idx, sign);
                        ptr += 1;
                    }
                }
            }

            if ptr == batch_size {
                batch_add2(batch_size, &mut a_bucks, &sched_bases, bases);
                Sched::clear(&mut sched_bases);

                ptr = 0;
            }
        }

        {
            batch_add2(ptr, &mut a_bucks, &sched_bases, bases);
            Sched::clear(&mut sched_bases);
            ptr = 0;
        }

        let mut running_sum = G1::identity();
        for (j_buck, a_buck) in j_bucks.into_iter().zip(a_bucks.into_iter()).rev() {
            running_sum += a_buck + j_buck;
            acc += &running_sum;
        }
    }

    // stat.debug();

    acc
}

pub fn multiexp_serial_w_booth(
    coeffs: &[Fr],
    bases: &[G1Affine],
    c: usize,
    batch_size: usize,
) -> G1 {
    let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();

    let segments = (256 / c) + 1;
    let mut acc = G1::identity();

    // let mut stat = Stat::default();

    for current_segment in (0..segments).rev() {
        for _ in 0..c {
            acc = acc.double();
        }

        let mut sched_bases = vec![Sched::default(); batch_size];
        let mut ptr = 0;

        let mut a_bucks = vec![G1Affine::identity(); 1 << (c - 1)];
        let mut j_bucks = vec![G1::identity(); 1 << (c - 1)];

        for (base_idx, coeff) in coeffs.iter().enumerate() {
            let buck_idx = get_booth_index(current_segment, c, coeff.as_ref());

            if buck_idx != 0 {
                let sign = buck_idx.is_positive();
                let buck_idx = buck_idx.unsigned_abs() as usize - 1;

                if Sched::contains(&sched_bases, buck_idx as usize) {
                    // stat.add(current_segment, 0, 0, 1);
                    if sign {
                        j_bucks[buck_idx] += bases[base_idx];
                    } else {
                        j_bucks[buck_idx] += bases[base_idx].neg();
                    }
                } else {
                    // stat.add(current_segment, 0, 1, 0);
                    if bool::from(a_bucks[buck_idx].is_identity()) {
                        if sign {
                            a_bucks[buck_idx] = bases[base_idx];
                        } else {
                            a_bucks[buck_idx] = bases[base_idx].neg();
                        }
                    } else {
                        sched_bases[ptr] = Sched::new(base_idx, buck_idx, sign);
                        ptr += 1;
                    }
                }
            }

            if ptr == batch_size {
                Sched::batch_add(&mut sched_bases, batch_size, &mut a_bucks, bases);
                ptr = 0;
            }
        }

        Sched::batch_add(&mut sched_bases, ptr, &mut a_bucks, bases);

        let mut running_sum = G1::identity();
        for (j_buck, a_buck) in j_bucks.into_iter().zip(a_bucks.into_iter()).rev() {
            running_sum += a_buck + j_buck;
            acc += &running_sum;
        }
    }

    // stat.debug();

    acc
}

pub struct MSM {
    sched_bases: Vec<Sched>,
    batch_size: usize,
}

impl MSM {
    pub fn init(batch_size: usize) -> Self {
        Self {
            sched_bases: vec![Sched::default(); batch_size],
            batch_size,
        }
    }
    pub fn multiexp_serial(&mut self, coeffs: &[Fr], bases: &[G1Affine], c: usize) -> G1 {
        let t0 = start_timer!(|| "to_repr");
        let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();
        end_timer!(t0);

        let segments = (256 / c) + 1;
        let mut acc = G1::identity();

        let mut ptr = 0;

        let mut a_bucks = vec![G1Affine::default(); 1 << (c - 1)];
        let mut j_bucks = vec![G1::default(); 1 << (c - 1)];

        for current_segment in (0..segments).rev() {
            for _ in 0..c {
                acc = acc.double();
            }

            for (base_idx, coeff) in coeffs.iter().enumerate() {
                let buck_idx = get_booth_index(current_segment, c, coeff.as_ref());

                if buck_idx != 0 {
                    let sign = buck_idx.is_positive();
                    let buck_idx = buck_idx.unsigned_abs() as usize - 1;

                    if Sched::contains(&self.sched_bases, buck_idx as usize) {
                        if sign {
                            j_bucks[buck_idx] += bases[base_idx];
                        } else {
                            j_bucks[buck_idx] -= bases[base_idx];
                        }
                    } else {
                        if bool::from(a_bucks[buck_idx].is_identity()) {
                            if sign {
                                a_bucks[buck_idx] = bases[base_idx];
                            } else {
                                a_bucks[buck_idx] = bases[base_idx].neg();
                            }
                        } else {
                            self.sched_bases[ptr] = Sched::new(base_idx, buck_idx, sign);
                            ptr += 1;
                        }
                    }
                }

                if ptr == self.batch_size {
                    batch_add2(self.batch_size, &mut a_bucks, &self.sched_bases, bases);
                    ptr = 0;
                }
            }

            {
                batch_add2(ptr, &mut a_bucks, &self.sched_bases, bases);
                ptr = 0;
            }

            let mut running_sum = G1::identity();
            for (j_buck, a_buck) in j_bucks.iter().zip(a_bucks.iter()).rev() {
                running_sum += a_buck + j_buck;
                acc += &running_sum;
            }
            a_bucks
                .iter_mut()
                .for_each(|a_buck| *a_buck = G1Affine::identity());
            j_bucks
                .iter_mut()
                .for_each(|j_buck| *j_buck = G1::identity());
        }

        acc
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

    fn affine(&self) -> G1Affine {
        match *self {
            Bucket::None => G1Affine::identity(),
            Bucket::Point(a) => a.to_affine(),
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
            Bucket::None => panic!("None"),
            Bucket::Point(a) => a.x,
        }
    }

    fn y(&self) -> Fq {
        match *self {
            Bucket::None => panic!("None"),
            Bucket::Point(a) => a.y,
        }
    }

    fn set_x(&mut self, x: &Fq) {
        match *self {
            Bucket::None => panic!("None"),
            Bucket::Point(ref mut a) => a.x = *x,
        }
    }

    // fn is_on_curve(&mut self) {
    //     match *self {
    //         Bucket::None => panic!("None"),
    //         Bucket::Point(ref mut a) => assert!(bool::from(a.is_on_curve())),
    //     }
    // }

    fn set_y(&mut self, y: &Fq) {
        match *self {
            Bucket::None => panic!("None"),
            Bucket::Point(ref mut a) => a.y = *y,
        }
    }
}

pub fn multiexp_serial_non_bucks(
    coeffs: &[Fr],
    bases: &[G1Affine],
    c: usize,
    batch_size: usize,
) -> G1 {
    let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();

    let segments = (256 / c) + 1;
    let mut acc = G1::identity();

    // let mut stat = Stat::default();

    for current_segment in (0..segments).rev() {
        for _ in 0..c {
            acc = acc.double();
        }

        let mut a_bucks = vec![Bucket::None; 1 << (c - 1)];
        let mut j_bucks = vec![Bucket::None; 1 << (c - 1)];

        let mut sched_bases = vec![Sched::default(); batch_size];
        let mut ptr = 0;

        for (base_idx, coeff) in coeffs.iter().enumerate() {
            let buck_idx = get_booth_index(current_segment, c, coeff.as_ref());

            if buck_idx != 0 {
                let sign = buck_idx.is_positive();
                let buck_idx = buck_idx.unsigned_abs() as usize - 1;

                if Sched::contains(&sched_bases, buck_idx as usize) {
                    // stat.add(current_segment, 0, 0, 1);
                    j_bucks[buck_idx].add_assign(&bases[base_idx], sign);
                } else {
                    // stat.add(current_segment, 0, 1, 0);
                    if !a_bucks[buck_idx].assign(&bases[base_idx], sign) {
                        sched_bases[ptr] = Sched::new(base_idx, buck_idx, sign);
                        ptr += 1;
                    }
                }
            }

            if ptr == batch_size {
                batch_add3(batch_size, &mut a_bucks[..], &sched_bases, bases);
                Sched::clear(&mut sched_bases);
                ptr = 0;
            }
        }

        batch_add3(ptr, &mut a_bucks[..], &sched_bases, bases);
        Sched::clear(&mut sched_bases);

        let mut running_sum = G1::identity();
        for (j_buck, a_buck) in j_bucks.iter().zip(a_bucks.iter()).rev() {
            running_sum += j_buck.add(a_buck);
            acc += &running_sum;
        }
    }
    acc
}

pub fn multiexp_serial_second_pass(
    coeffs: &[Fr],
    bases: &[G1Affine],
    c: usize,
    batch_size: usize,
) -> G1 {
    let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();

    let segments = (256 / c) + 1;
    let mut acc = G1::identity();

    // let mut stat = Stat::default();

    for current_segment in (0..segments).rev() {
        for _ in 0..c {
            acc = acc.double();
        }

        let mut a_bucks = vec![Bucket::None; 1 << (c - 1)];
        let mut j_bucks = vec![Bucket::None; 1 << (c - 1)];

        let mut sched_bases0 = vec![Sched::default(); batch_size];
        let mut sched_bases1 = vec![Sched::default(); batch_size];
        let mut ptr0 = 0;
        let mut ptr1 = 0;

        for (base_idx, coeff) in coeffs.iter().enumerate() {
            let buck_idx = get_booth_index(current_segment, c, coeff.as_ref());

            if buck_idx != 0 {
                let sign = buck_idx.is_positive();
                let buck_idx = buck_idx.unsigned_abs() as usize - 1;

                if Sched::contains(&sched_bases0, buck_idx as usize) {
                    if Sched::contains(&sched_bases1, buck_idx as usize) {
                        // stat.add(current_segment, 0, 0, 1);
                        j_bucks[buck_idx].add_assign(&bases[base_idx], sign);
                    } else {
                        // stat.add(current_segment, 0, 1, 0);
                        if !a_bucks[buck_idx].assign(&bases[base_idx], sign) {
                            sched_bases1[ptr1] = Sched::new(base_idx, buck_idx, sign);
                            ptr1 += 1;
                        }
                    }
                } else {
                    // stat.add(current_segment, 0, 1, 0);
                    if !a_bucks[buck_idx].assign(&bases[base_idx], sign) {
                        sched_bases0[ptr0] = Sched::new(base_idx, buck_idx, sign);
                        ptr0 += 1;
                    }
                }
            }

            if ptr0 == batch_size {
                batch_add3(batch_size, &mut a_bucks[..], &sched_bases0, bases);
                Sched::clear(&mut sched_bases0);
                ptr0 = 0;
            }
            if ptr1 == batch_size {
                batch_add3(batch_size, &mut a_bucks[..], &sched_bases1, bases);
                Sched::clear(&mut sched_bases1);
                ptr1 = 0;
            }
        }

        batch_add3(ptr0, &mut a_bucks[..], &sched_bases0, bases);
        Sched::clear(&mut sched_bases0);

        batch_add3(ptr1, &mut a_bucks[..], &sched_bases1, bases);
        Sched::clear(&mut sched_bases1);

        let mut running_sum = G1::identity();
        for (j_buck, a_buck) in j_bucks.iter().zip(a_bucks.iter()).rev() {
            running_sum += j_buck.add(a_buck);
            acc += &running_sum;
        }
    }
    // stat.debug();
    acc
}

// pub struct MSM2 {
//     // sched_bases: Vec<Sched>,
//     // batch_size: usize,
//     // a_bucks: Vec<Bucket<G1Affine>>,
//     // j_bucks: Vec<Bucket<G1>>,
// }

// impl MSM2 {
//     pub fn init(batch_size: usize, c: usize) -> Self {
//         Self {
//             // sched_bases: vec![Sched::default(); batch_size],
//             // batch_size,
//             // a_bucks: vec![Bucket::None; 1 << (c - 1)],
//             // j_bucks: vec![Bucket::None; 1 << (c - 1)],
//         }
//     }
//     pub fn multiexp_serial(&mut self, coeffs: &[Fr], bases: &[G1Affine], c: usize) -> G1 {
//         let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();

//         let segments = (256 / c) + 1;
//         let mut acc = G1::identity();

//         // let mut stat = Stat::default();

//         for current_segment in (0..segments).rev() {
//             // println!("XXX {}", current_segment);
//             for _ in 0..c {
//                 acc = acc.double();
//             }

//             let mut a_bucks = vec![Bucket::None; 1 << (c - 1)];
//             let mut j_bucks = vec![Bucket::None; 1 << (c - 1)];

//             let mut sched_bases = vec![Sched::default(); 64];
//             let mut ptr = 0;

//             for (base_idx, coeff) in coeffs.iter().enumerate() {
//                 let buck_idx = get_booth_index(current_segment, c, coeff.as_ref());

//                 if buck_idx != 0 {
//                     let sign = buck_idx.is_positive();
//                     let buck_idx = buck_idx.unsigned_abs() as usize - 1;

//                     if Sched::contains(&sched_bases, buck_idx as usize) {
//                         // stat.add(current_segment, 0, 0, 1);
//                         j_bucks[buck_idx].add_assign(&bases[base_idx], sign);
//                     } else {
//                         // stat.add(current_segment, 0, 1, 0);
//                         if !a_bucks[buck_idx].assign(&bases[base_idx], sign) {
//                             sched_bases[ptr] = Sched::new(base_idx, buck_idx, sign);
//                             ptr += 1;
//                         }
//                     }
//                 }

//                 if ptr == 64 {
//                     batch_add3(64, &mut a_bucks[..], &sched_bases, bases);
//                     Sched::clear(&mut sched_bases);
//                     ptr = 0;
//                 }
//             }

//             batch_add3(ptr, &mut a_bucks[..], &sched_bases, bases);
//             Sched::clear(&mut sched_bases);

//             let mut running_sum = G1::identity();
//             for (j_buck, a_buck) in j_bucks.iter().zip(a_bucks.iter()).rev() {
//                 running_sum += j_buck.add(a_buck);
//                 acc += &running_sum;
//             }

//             // a_bucks.iter_mut().for_each(|a_buck| *a_buck = Bucket::None);
//             // j_bucks.iter_mut().for_each(|j_buck| *j_buck = Bucket::None);
//         }

//         // stat.debug();

//         acc
//     }
// }

// pub fn multiexp_serial4(coeffs: &[Fr], bases: &[G1Affine], c: usize, batch_size: usize) -> G1 {
//     let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();

//     let segments = (256 / c) + 1;
//     let mut acc = G1::identity();

//     let mut sched_points = vec![Sched::default(); batch_size];
//     let mut ptr = 0;

//     #[derive(Debug, Clone)]
//     enum Bucket {
//         None,
//         Point(Option<G1Affine>, Option<G1>),
//     }

//     impl Bucket {
//         fn is_identity(&self) -> bool {
//             match self {
//                 Bucket::None => true,
//                 _ => false,
//             }
//         }
//     }

//     for current_segment in (0..segments).rev() {
//         for _ in 0..c {
//             acc = acc.double();
//         }

//         let mut bucks = vec![Bucket::None; 1 << (c - 1)];
//         let mut a_bucks = vec![G1Affine::identity(); 1 << (c - 1)];
//         let mut j_bucks = vec![G1::identity(); 1 << (c - 1)];

//         for (base_idx, coeff) in coeffs.iter().enumerate() {
//             let buck_idx = get_booth_index(current_segment, c, coeff.as_ref());

//             if buck_idx != 0 {
//                 let sign = buck_idx.is_positive();
//                 let buck_idx = buck_idx.unsigned_abs() as usize - 1;

//                 // match &mut bucks[buck_idx] {
//                 //     Bucket::None => {
//                 //         bucks[buck_idx] = Bucket::Point(Some(bases[base_idx]), None);
//                 //     }
//                 //     Bucket::Point(affine, jack) => {
//                 //         //
//                 //         affine.map(|_| {
//                 //             *j_buck += if sign {
//                 //                 bases[base_idx]
//                 //             } else {
//                 //                 bases[base_idx].neg()
//                 //             }
//                 //         })
//                 //     }
//                 //     Bucket::Sched(a_buck) => {
//                 //         //
//                 //     }
//                 // }

//                 if Sched::contains(&sched_points, buck_idx as usize) {
//                     if sign {
//                         j_bucks[buck_idx] += bases[base_idx];
//                     } else {
//                         j_bucks[buck_idx] -= bases[base_idx];
//                     }
//                 } else {
//                     if bool::from(a_bucks[buck_idx].is_identity()) {
//                         if sign {
//                             a_bucks[buck_idx] = bases[base_idx];
//                         } else {
//                             a_bucks[buck_idx] = bases[base_idx].neg();
//                         }
//                     } else {
//                         sched_points[ptr] = Sched::new(base_idx, buck_idx, sign);
//                         ptr += 1;
//                     }
//                 }
//             }

//             if ptr == batch_size {
//                 // let t0 = start_timer!(|| "batch_add2");
//                 batch_add2(batch_size, &mut a_bucks, &sched_points, bases);
//                 // end_timer!(t0);
//                 ptr = 0;
//             }
//         }

//         {
//             batch_add2(ptr, &mut a_bucks, &sched_points, bases);
//             ptr = 0;
//         }

//         let mut running_sum = G1::identity();
//         for (j_buck, a_buck) in j_bucks.into_iter().zip(a_bucks.into_iter()).rev() {
//             running_sum += a_buck + j_buck;
//             acc += &running_sum;
//         }
//     }

//     acc
// }

// sched_bases
//     .iter()
//     .position(|(sched_idx, _)| *sched_idx == buck_idx)
//     .map_or_else(
//         || {
//             if bool::from(a_bucks[buck_idx].is_identity()) {
//                 a_bucks[buck_idx] = bases[base_idx];
//             } else {
//                 sched_bases[ptr] = (buck_idx, base_idx);
//                 ptr += 1;
//             }
//         },
//         |_| {
//             j_bucks[buck_idx] += bases[base_idx];
//         },
//     );
// {
//     let mut t = vec![Fq::ZERO; batch_size];
//     let mut acc = Fq::ONE;

//     for ((buck_idx, base_idx), t) in sched_bucks.iter().zip(t.iter_mut()) {
//         *t = acc;
//         acc *= a_bucks[*buck_idx].x - bases[*base_idx].x;
//     }

//     acc = acc.invert().unwrap();

//     for ((buck_idx, base_idx), t) in sched_bucks.iter().zip(t.iter()).rev() {
//         let inv = acc * t;
//         let base = bases[*base_idx];
//         acc = acc * (a_bucks[*buck_idx].x - base.x);

//         let lambda = (a_bucks[*buck_idx].y - base.y) * inv;
//         let x = lambda.square() - (a_bucks[*buck_idx].x + base.x);
//         a_bucks[*buck_idx].y = (lambda * (base.x - x)) - base.y;
//         a_bucks[*buck_idx].x = x;
//     }
// }

// println!(
//     "segment: {}, n_jac: {}, n_aff: {}",
//     current_segment, stat.n_jac, stat.n_aff
// );

// j_bucks
//     .iter_mut()
//     .zip(a_bucks.iter())
//     .zip(init_buckets.iter())
//     .for_each(|((j_buck, a_buck), init)| {
//         *j_buck += a_buck - init;
//     });
