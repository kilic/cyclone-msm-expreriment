use std::collections::BTreeMap;

pub mod cyclone;
pub mod cyclone_dev;
pub mod cyclone_par_by_window;
pub mod pr29;
#[cfg(test)]
mod test;
pub mod upstream;
pub mod utils;

#[derive(Debug, Default)]
pub(crate) struct Stat(BTreeMap<usize, (usize, usize, usize)>);

impl Stat {
    pub(crate) fn add(&mut self, seg: usize, n_zero: usize, n_aff: usize, n_jac: usize) {
        self.0
            .entry(seg)
            .and_modify(|(zero, aff, jac)| {
                *aff += n_aff;
                *jac += n_jac;
                *zero += n_zero;
            })
            .or_insert((n_zero, n_aff, n_aff));
    }

    pub(crate) fn debug(&self) {
        let mut tot_aff = 0;
        let mut tot_jac = 0;
        let mut tot_zero = 0;
        for (_, (n_zero, n_aff, n_jac)) in self.0.iter() {
            tot_aff += n_aff;
            tot_jac += n_jac;
            tot_zero += n_zero;
        }
        println!(
            "tot_aff: {}, tot_aff: {}, tot_jac: {}",
            tot_zero, tot_aff, tot_jac
        )
    }
}
