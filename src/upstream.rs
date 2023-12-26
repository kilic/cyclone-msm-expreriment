pub mod zcash {
    use crate::utils::BucketIndex;
    use ff::PrimeField;
    use group::Group;
    use halo2curves::CurveAffine;
    use rayon::{current_num_threads, scope};

    pub fn msm_serial<C: CurveAffine>(coeffs: &[C::Scalar], bases: &[C]) -> C::Curve {
        let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();

        let c = if bases.len() < 4 {
            1
        } else if bases.len() < 32 {
            3
        } else {
            (f64::from(bases.len() as u32)).ln().ceil() as usize
        };

        let segments = (256 / c) + 1;

        let mut acc = C::Curve::identity();

        for current_segment in (0..segments).rev() {
            for _ in 0..c {
                acc = acc.double();
            }

            #[derive(Clone, Copy)]
            enum Bucket<C: CurveAffine> {
                None,
                Affine(C),
                Projective(C::Curve),
            }

            impl<C: CurveAffine> Bucket<C> {
                fn add_assign(&mut self, other: &C) {
                    *self = match *self {
                        Bucket::None => Bucket::Affine(*other),
                        Bucket::Affine(a) => Bucket::Projective(a + *other),
                        Bucket::Projective(mut a) => {
                            a += *other;
                            Bucket::Projective(a)
                        }
                    }
                }

                fn add(self, mut other: C::Curve) -> C::Curve {
                    match self {
                        Bucket::None => other,
                        Bucket::Affine(a) => {
                            other += a;
                            other
                        }
                        Bucket::Projective(a) => other + a,
                    }
                }
            }

            let mut buckets: Vec<Bucket<_>> = vec![Bucket::None; (1 << c) - 1];

            for (coeff, base) in coeffs.iter().zip(bases.iter()) {
                let coeff = crate::utils::Slice::get(current_segment, c, coeff.as_ref());
                if coeff != 0 {
                    buckets[coeff - 1].add_assign(base);
                }
            }

            // Summation by parts
            // e.g. 3a + 2b + 1c = a +
            //                    (a) + b +
            //                    ((a) + b) + c
            let mut running_sum = C::Curve::identity();
            for exp in buckets.into_iter().rev() {
                running_sum = exp.add(running_sum);
                acc += &running_sum;
            }
        }
        acc
    }

    pub fn msm_par<C: CurveAffine>(coeffs: &[C::Scalar], bases: &[C]) -> C::Curve {
        assert_eq!(coeffs.len(), bases.len());

        let num_threads = current_num_threads();
        if coeffs.len() > num_threads {
            let chunk = coeffs.len() / num_threads;
            let num_chunks = coeffs.chunks(chunk).len();
            let mut results = vec![C::Curve::identity(); num_chunks];
            scope(|scope| {
                let chunk = coeffs.len() / num_threads;

                for ((coeffs, bases), acc) in coeffs
                    .chunks(chunk)
                    .zip(bases.chunks(chunk))
                    .zip(results.iter_mut())
                {
                    scope.spawn(move |_| {
                        *acc = msm_serial(coeffs, bases);
                    });
                }
            });
            results.iter().fold(C::Curve::identity(), |a, b| a + b)
        } else {
            msm_serial(coeffs, bases)
        }
    }
}

pub mod pse {

    use ff::PrimeField;
    use group::Group;
    use halo2curves::CurveAffine;
    use rayon::{current_num_threads, scope};

    use crate::utils::BucketIndex;

    pub fn msm_serial<C: CurveAffine>(coeffs: &[C::Scalar], bases: &[C]) -> C::Curve {
        let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();

        let c = if bases.len() < 4 {
            1
        } else if bases.len() < 32 {
            3
        } else {
            (f64::from(bases.len() as u32)).ln().ceil() as usize
        };

        let number_of_windows = C::Scalar::NUM_BITS as usize / c + 1;

        let mut acc = C::Curve::identity();

        for current_window in (0..number_of_windows).rev() {
            for _ in 0..c {
                acc = acc.double();
            }

            #[derive(Clone, Copy)]
            enum Bucket<C: CurveAffine> {
                None,
                Affine(C),
                Projective(C::Curve),
            }

            impl<C: CurveAffine> Bucket<C> {
                fn add_assign(&mut self, other: &C) {
                    *self = match *self {
                        Bucket::None => Bucket::Affine(*other),
                        Bucket::Affine(a) => Bucket::Projective(a + *other),
                        Bucket::Projective(mut a) => {
                            a += *other;
                            Bucket::Projective(a)
                        }
                    }
                }

                fn add(self, mut other: C::Curve) -> C::Curve {
                    match self {
                        Bucket::None => other,
                        Bucket::Affine(a) => {
                            other += a;
                            other
                        }
                        Bucket::Projective(a) => other + a,
                    }
                }
            }

            let mut buckets: Vec<Bucket<C>> = vec![Bucket::None; 1 << (c - 1)];

            for (coeff, base) in coeffs.iter().zip(bases.iter()) {
                let coeff = crate::utils::Booth::get(current_window, c, coeff.as_ref());
                if coeff.is_positive() {
                    buckets[coeff as usize - 1].add_assign(base);
                }
                if coeff.is_negative() {
                    buckets[coeff.abs() as usize - 1].add_assign(&base.neg());
                }
            }

            // Summation by parts
            // e.g. 3a + 2b + 1c = a +
            //                    (a) + b +
            //                    ((a) + b) + c
            let mut running_sum = C::Curve::identity();
            for exp in buckets.into_iter().rev() {
                running_sum = exp.add(running_sum);
                acc += &running_sum;
            }
        }

        acc
    }

    pub fn msm_par<C: CurveAffine>(coeffs: &[C::Scalar], bases: &[C]) -> C::Curve {
        assert_eq!(coeffs.len(), bases.len());

        let num_threads = current_num_threads();
        if coeffs.len() > num_threads {
            let chunk = coeffs.len() / num_threads;
            let num_chunks = coeffs.chunks(chunk).len();
            let mut results = vec![C::Curve::identity(); num_chunks];
            scope(|scope| {
                let chunk = coeffs.len() / num_threads;

                for ((coeffs, bases), acc) in coeffs
                    .chunks(chunk)
                    .zip(bases.chunks(chunk))
                    .zip(results.iter_mut())
                {
                    scope.spawn(move |_| {
                        *acc = msm_serial(coeffs, bases);
                    });
                }
            });
            results.iter().fold(C::Curve::identity(), |a, b| a + b)
        } else {
            msm_serial(coeffs, bases)
        }
    }
}
