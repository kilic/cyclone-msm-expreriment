use std::env::consts::OS;

use ark_std::{end_timer, start_timer};
use ff::{Field, PrimeField};
use group::{prime::PrimeCurveAffine, Curve, Group};
use halo2curves::{
    bn256::Fr,
    bn256::{G1Affine, G1},
};
use rand_core::OsRng;

fn get_data(k: usize) -> (Vec<G1Affine>, Vec<Fr>) {
    // TODO: write to a file
    let points = (0..1 << k).map(|_| G1::random(OsRng)).collect::<Vec<_>>();
    let mut affine_points = vec![G1Affine::identity(); 1 << k];
    G1::batch_normalize(&points[..], &mut affine_points[..]);
    let scalars = (0..1 << k).map(|_| Fr::random(OsRng)).collect::<Vec<_>>();
    (affine_points, scalars)
}

#[test]
fn test_cyclone() {
    let (min_k, max_k) = (18, 19);
    let (points, scalars) = get_data(max_k);

    for k in min_k..=max_k {
        println!("k: {}", k);

        let points = &points[..1 << k];
        let scalars = &scalars[..1 << k];

        let t0 = start_timer!(|| format!("zcash serial"));
        let e0 = crate::upstream::zcash::msm_serial(scalars, points);
        end_timer!(t0);

        let t0 = start_timer!(|| format!("* zcash par"));
        let e1 = crate::upstream::zcash::msm_par(scalars, points);
        end_timer!(t0);
        assert_eq!(e0, e1);

        let t0 = start_timer!(|| format!("pse serial"));
        let e1 = crate::upstream::pse::msm_serial(scalars, points);
        end_timer!(t0);
        assert_eq!(e0, e1);

        let t0 = start_timer!(|| format!("* pse par"));
        let e1 = crate::upstream::pse::msm_par(scalars, points);
        end_timer!(t0);
        assert_eq!(e0, e1);

        let t0 = start_timer!(|| format!("* pr29"));
        let e1 = crate::pr29::MSM::best(scalars, points);
        end_timer!(t0);
        assert_eq!(e0, e1);

        let (min_c, max_c) = (10, 14);
        let batch_size = 64;
        for c in min_c..=max_c {
            println!("c: {}", c);
            // println!("mamy says batch size {} ", 4 * c * c - 16 * c - 128);

            let t0 = start_timer!(|| format!("cyclone serial {} {}", c, batch_size));
            let e1 = crate::cyclone::msm_serial(scalars, points, c, batch_size);
            end_timer!(t0);
            assert_eq!(e0, e1);

            let t0 = start_timer!(|| format!("cyclone par split {} {}", c, batch_size));
            let e1 = crate::cyclone::msm_par_split(scalars, points, c, batch_size);
            end_timer!(t0);
            assert_eq!(e0, e1);

            let t0 = start_timer!(|| format!("* cyclone par window {} {}", c, batch_size));
            let e1 = crate::cyclone::msm_par_window(scalars, points, c, batch_size);
            end_timer!(t0);
            assert_eq!(e0, e1);
        }
        println!("---");
    }
}
