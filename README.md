Experimentations of MSM with batch addition and scheduler idea from [CycloneMSM](https://eprint.iacr.org/2022/1396.pdf) using halo2curves as backend

Previous efforts:
  * [#40](https://github.com/privacy-scaling-explorations/halo2/pull/40)
  * [#29](https://github.com/privacy-scaling-explorations/halo2curves/pull/29)


This should be faster than those and hopefully close to gnark-crypto level and much friendly to review eventually.

### Notes 
* for now only with greedy scheduler
* parallelizing experiments:
  * [ ] chunk into smaller MSMs as usual
  * [ ] chunk by buckets as in [constantine](https://github.com/mratsim/constantine/blob/master/constantine/math/elliptic/ec_multi_scalar_mul_scheduler.nim)
  * [x] chunk by window as in [gnark-crypto](https://github.com/Consensys/gnark-crypto/blob/master/ecc/bn254/multiexp_affine.go) [ashWhiteHat-PR](https://github.com/zcash/halo2/pull/796)
* [zkcrypto](https://github.com/zkcrypto/group) RFC: mutable access to `CurveAffine` coordinates without `is_on_curve` check. Due to no such support it is bn256 only

