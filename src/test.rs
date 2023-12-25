use std::env::consts::OS;

use ark_std::{end_timer, start_timer};
use ff::{Field, PrimeField};
use group::{prime::PrimeCurveAffine, Curve, Group};
use halo2curves::{
    bn256::Fr,
    bn256::{G1Affine, G1},
};
use rand_core::OsRng;

use crate::{
    pr29,
    utils::{get_booth_index, get_bucket_index},
};

fn get_data(k: usize) -> (Vec<G1Affine>, Vec<Fr>) {
    let points = (0..1 << k).map(|_| G1::random(OsRng)).collect::<Vec<_>>();
    let mut affine_points = vec![G1Affine::identity(); 1 << k];
    G1::batch_normalize(&points[..], &mut affine_points[..]);
    let scalars = (0..1 << k).map(|_| Fr::random(OsRng)).collect::<Vec<_>>();
    (affine_points, scalars)
}

#[test]
fn test_inv() {}

#[test]
fn test_cyclone() {
    let (min_k, max_k) = (20, 20);
    let (points, scalars) = get_data(max_k);

    for k in min_k..=max_k {
        println!("k: {}", k);

        let points = &points[..1 << k];
        let scalars = &scalars[..1 << k];

        let t0 = start_timer!(|| format!("zcash serial"));
        let e0 = crate::upstream::zcash::multiexp_serial(scalars, points);
        end_timer!(t0);

        let t0 = start_timer!(|| format!("zcash parr"));
        let e1 = crate::upstream::zcash::best_multiexp(scalars, points);
        end_timer!(t0);
        assert_eq!(e0, e1);

        let t0 = start_timer!(|| format!("pse serial"));
        let e1 = crate::upstream::zcash::multiexp_serial(scalars, points);
        end_timer!(t0);
        assert_eq!(e0, e1);

        let t0 = start_timer!(|| format!("pse parr"));
        let e1 = crate::upstream::zcash::best_multiexp(scalars, points);
        end_timer!(t0);
        assert_eq!(e0, e1);

        // let t0 = start_timer!(|| format!("zcash parr"));
        // let e0 = crate::zcash::best_multiexp(scalars, points);
        // end_timer!(t0);
        // println!("e0: {:?}", e0.to_affine());

        let (min_c, max_c) = (10, 16);
        let batch_size = 64;
        for c in min_c..=max_c {
            println!("c: {}", c);
            // println!("mamy says {} ", 4 * c * c - 16 * c - 128);

            let t0 = start_timer!(|| format!("cyclone {} {}", c, batch_size));
            let e1 = crate::cyclone::multiexp_serial(scalars, points, c, batch_size);
            end_timer!(t0);
            assert_eq!(e0, e1);

            // let t0 = start_timer!(|| format!("cyclone2 {} {}", c, batch_size));
            // let e1 = crate::cyclone::multiexp_serial2(scalars, points, c, batch_size);
            // end_timer!(t0);
            // assert_eq!(e0, e1);

            // let t0 = start_timer!(|| format!("cyclone332 {} {}", c, batch_size));
            // let e1 = crate::cyclone::multiexp_serial3(scalars, points, c, 32);
            // end_timer!(t0);
            // assert_eq!(e0, e1);

            // let t0 = start_timer!(|| format!("cyclone3 {} {}", c, batch_size));
            // let e1 = crate::cyclone::multiexp_serial3(scalars, points, c, batch_size);
            // end_timer!(t0);
            // assert_eq!(e0, e1);

            // let t0 = start_timer!(|| format!("cyclone4 {} {}", c, batch_size));
            // let e1 = crate::cyclone::multiexp_serial4(scalars, points, c, 64);
            // end_timer!(t0);
            // assert_eq!(e0, e1);

            let t0 = start_timer!(|| format!("cyclone-w-b {} {}", c, batch_size));
            let e1 = crate::cyclone::multiexp_serial_w_booth(scalars, points, c, batch_size);
            end_timer!(t0);
            assert_eq!(e0, e1);

            let t0 = start_timer!(|| format!("cyclone-nob {} {}", c, batch_size));
            let e1 = crate::cyclone::multiexp_serial_non_bucks(scalars, points, c, batch_size);
            end_timer!(t0);
            assert_eq!(e0, e1);

            let t0 = start_timer!(|| format!("cyclone-new {} {}", c, batch_size));
            let e1 = crate::cyclone_dev::multiexp_serial(scalars, points, c, batch_size);
            end_timer!(t0);
            assert_eq!(e0, e1);

            let t0 = start_timer!(|| format!("cyclone-par {} {}", c, batch_size));
            let e1 = crate::cyclone_par_by_window::multiexp_par(scalars, points, c, batch_size);
            end_timer!(t0);
            assert_eq!(e0, e1);

            // let t0 = start_timer!(|| format!("cyclone-secp {} {}", c, batch_size));
            // let e1 = crate::cyclone::multiexp_serial_second_pass(scalars, points, c, batch_size);
            // end_timer!(t0);
            // assert_eq!(e0, e1);

            // let mut msm = crate::cyclone::MSM3::init(batch_size, c);
            // let t0 = start_timer!(|| format!("cyclone-strct {} {}", c, batch_size));
            // let e1 = msm.multiexp_serial(scalars, points, c);
            // end_timer!(t0);
            // assert_eq!(e0, e1);

            // let t0 = start_timer!(|| format!("cyclone4 {} {}", c, 128));
            // let e1 = crate::cyclone::multiexp_serial4(scalars, points, c, 128);
            // end_timer!(t0);
            // assert_eq!(e0, e1);

            // let t0 = start_timer!(|| format!("cyclone4 {} {}", c, 256));
            // let e1 = crate::cyclone::multiexp_serial4(scalars, points, c, 256);
            // end_timer!(t0);
            // assert_eq!(e0, e1);

            // let mut msm = crate::cyclone::MSM::init(batch_size);
            // let t0 = start_timer!(|| format!("cyclones {} {}", c, batch_size));
            // let e1 = msm.multiexp_serial3(scalars, points, c);
            // end_timer!(t0);
            // assert_eq!(e0, e1);

            // let t0 = start_timer!(|| format!("cyclone128 {} {}", c, batch_size));
            // let e1 = crate::cyclone::multiexp_serial3(scalars, points, c, 128);
            // end_timer!(t0);
            // assert_eq!(e0, e1);

            // for batch_k in 5..=12 {

            //     // assert_eq!(e0, e1);
            // }
            // println!("e1: {:?}", e1.to_affine());
        }

        // let t0 = start_timer!(|| format!("kilic"));

        // let mut e2 = G1::identity();
        // kilic29::MSM::evaluate_with(scalars, points, &mut e2, None);
        // end_timer!(t0);

        // println!("{}", e2 == e0);
    }
}
