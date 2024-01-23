MSM Experimentations with batch addition and scheduler idea from [CycloneMSM](https://eprint.iacr.org/2022/1396.pdf) combined with suggestions at [halo2/#187](https://github.com/privacy-scaling-explorations/halo2/issues/187)

Previous efforts that targets `halo2curves`:

* with batch addition and index ordering as in barretenberg approach:
  * [pse/halo2/#40](https://github.com/privacy-scaling-explorations/halo2/pull/40)
  * [pse/halo2curves/#29](https://github.com/privacy-scaling-explorations/halo2curves/pull/29)
* with chunk by window approach
  * [ashWhiteHat-PR](https://github.com/zcash/halo2/pull/796)


Goal: This should be faster than those even in high `k` param and hopefully close to gnark-crypto level and eventually becomes much friendly for review.

### Notes 
* For now only greedy scheduler approach is implemented. Not only that we are leaving bunch of surfaces to refine.
* Parallelizing experiments:
  * [x] chunk into smaller MSMs as usual
  * [ ] chunk by buckets as in [constantine](https://github.com/mratsim/constantine/blob/master/constantine/math/elliptic/ec_multi_scalar_mul_scheduler.nim)
  * [x] chunk by window as in [gnark-crypto](https://github.com/Consensys/gnark-crypto/blob/master/ecc/bn254/multiexp_affine.go) [ashWhiteHat-PR](https://github.com/zcash/halo2/pull/796)
* [zkcrypto](https://github.com/zkcrypto/group) RFC: mutable access to coordinates of `CurveAffine` without `is_on_curve` check. Without such support we can't go generic.
* Without this speed up [Fast modular inverse](https://github.com/privacy-scaling-explorations/halo2curves/pull/83) this work would make no sense.


### TODO

* [ ] Impl delayed scheduler
* [ ] Auto selection of window
* [ ] Auto selection of batch size
* [ ] Chunk by bucket challenge

### Test

Run the test with timer:

```
cargo test cyclone --release -- --nocapture 
```

#### Local results:

| k   | #29   | current | this  | this w/ coods api |
| --- | ----- | -----   | ----- | ----- |
| 16  | 54ms  | 98ms    | 48ms  | 54ms  |
| 17  | 95ms  | 170ms   | 85ms  | 103ms |
| 18  | 175ms | 311ms   | 168ms | 185ms |
| 19  | 380ms | 577ms   | 287ms | 338ms |
| 20  | 760ms | 1.07s   | 567ms | 639ms |
| 21  | 1.58s | 1.98s   | 1.08s | 1.26s |
| 22  | 2.80s | 3.80s   | 2.14s | 2.49s |

run on M1 machine
