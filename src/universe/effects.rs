pub(crate) mod general_relativity;
pub(crate) mod magnetism;
pub(crate) mod tides;
pub(crate) mod wind;
pub use general_relativity::GeneralRelativityModel;
pub use magnetism::MagneticModel;
pub use tides::{Kaula, SpectrumFile};
pub use wind::WindModel;
