import logging
def main(input_mtz, input_pdb, intensity_column, sigma_column):
    import gemmi
    try:
        log = logging.getLogger(__name__)
        logging.basicConfig(level=logging.INFO)
        mtz = gemmi.read_mtz_file(input_mtz)
        mtz.update_reso()

        f_obs = mtz.get_value_sigma(intensity_column, sigma_column)
        struct = gemmi.read_structure(input_pdb)
        d_min = mtz.resolution_high() - 1e-9
        
        density_calc = gemmi.DensityCalculatorE()
        density_calc.d_min = d_min
        density_calc.set_refmac_compatible_blur(struct[0])
        density_calc.grid.setup_from(struct)
        density_calc.put_model_density_on_grid(struct[0])
        model_grid = gemmi.transform_map_to_f_phi(density_calc.grid, half_l=True)
        f_cryst = model_grid.prepare_asu_data(dmin=d_min, unblur=density_calc.blur)

        scaler = gemmi.Scaling(struct.cell, struct.find_spacegroup())
        scaler.use_solvent = False
        scaler.prepare_points(f_cryst, f_obs)
        scaler.fit_isotropic_b_approximately()
        scaler.fit_parameters()
        
        r_factor = round(scaler.calculate_r_factor(),5)
        log.info(f"Calculated scaled R-factor for {input_pdb}: {r_factor}")
        return r_factor
    except:
        log.warning(f"Failed calculation of R-factor for conformer {input_pdb}. Continuing.")
        return 100
        