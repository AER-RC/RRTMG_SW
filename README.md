# RRTMG_SW: Shortwave Radiative Transfer Model for GCMs

---
**Contents**

1. [Introduction](#intro)
2. [Releases](#releases)
3. [Column Version](#column)
4. [GCM Version](#gcm)
5. [Contact](#contact)
6. [References](#ref)

## Introduction <a name="intro"></a>

This package contains the source code and sample makefiles necessary to run the latest version of RRTMG\_SW, a correlated k-distribution shortwave radiative transfer model developed at AER for application to GCMs. This version of RRTMG\_SW utilizes a two-stream radiative transfer method as implemented at ECMWF. This code has also been modified to utilize updated FORTRAN coding features. Two modes of operation are possible: 

1) RRTMG\_SW can be run as a [column model](https://github.com/AER-RC/RRTMG_SW#rrtmg_sw--column-version) using the [input files](https://github.com/AER-RC/RRTMG_SW#input-data) and [source modules](https://github.com/AER-RC/RRTMG_SW#source-code) described in the column version section, or 
2) it can be implemented as a [subroutine into an atmospheric general circulation model or single column model](https://github.com/AER-RC/RRTMG_SW#rrtmg_sw--gcm-version).

The version of RRTMG\_SW provided here utilizes a reduced complement of 112 *g*-points, which is half of the 224 *g*-points used in the standard RRTM\_SW, and a two-stream method for radiative transfer. Additional minor changes have been made to enhance computational performance. Total fluxes are accurate to within 1-2 W m<sup>-2</sup> relative to the standard RRTM\_SW (using DISORT) in clear sky and in the presence of aerosols and within 6 W m<sup>-2</sup> in overcast sky. RRTM\_SW with DISORT is itself accurate to within 2 W m<sup>-2</sup> of the data-validated multiple scattering model, CHARTS. Required absorption coefficient input data can be read in either from data stored within the code or from a netCDF file as selected in the makefile. 

This model can also utilize McICA, the Monte-Carlo Independent Column Approximation, to represent sub-grid scale cloud variability such as cloud fraction and cloud overlap. If the McICA option is selected to model a cloudy profile in column mode, then the model will run stochastically, and the output fluxes and heating rates will be an average over 200 samples. In GCM mode, the code will calcualte a single column per profile, and the statistical basis is provided by the spatial and temporal dimensions of the 3-D calculations. Several cloud overlap methods are available for partial cloudiness including maximum-random, exponential, and exponential-random. Without McICA, RRTMG\_SW is limited to clear sky or overcast cloud conditions.

For more information on the model, see the [Wiki Description Page](https://github.com/AER-RC/RRTMG_SW/wiki/Description).

## Current Release <a name="releases"></a>

[Version 5.0 is the latest version of the model](https://github.com/AER-RC/RRTMG_SW/releases/tag/v5.0)

Releases before Version 5.0 are not publicly available.

## RRTMG_SW : Column Version <a name="column"></a>

### DOCUMENTATION
The following text files (in the `doc` directory), along with this `README` provide information on release updates and on using and running RRTMG\_SW:

| File Name | Description |
| :--- | :--- |
| `release_notes.txt` | Code archive update information |
| `rrtmg_sw_instructions.txt` | Input instructions for files `INPUT_RRTM`, `IN_CLD_RRTM` and `IN_AER_RRTM` |

### SOURCE CODE
The following source files (in the `src` directory) must be used to run RRTMG\_SW in stand-alone mode as a column model (the utility files are stored separately in the `aer_rt_utils` directory):

| File Name | Description |
| :--- | :--- |
| `rrtmg_sw.1col.f90` | Main module
| `rrtmg_sw_cldprop.f90` | Calculation of cloud optical properties
| `rrtmg_sw_cldprmc.f90` | Calculation of cloud optical properties (McICA)
| `rrtmg_sw_init.f90` | RRTMG_SW initialization routine; reduces *g*-intervals from 224 to 112
| `rrtmg_sw_k_g.f90` | Absorption coefficient data file
| `rrtmg_sw_read_nc.f90` | Optional absorption coefficient data netCDF input
| `rrtmg_sw_reftra.f90` | Calculation of two-stream reflectivities and transmissivities
| `rrtmg_sw_setcoef.f90` | Set up routine
| `rrtmg_sw_spcvrt.f90` | Top subroutine for two-stream model
| `rrtmg_sw_spcvmc.f90` | Top subroutine for two-stream model (McICA)
| `rrtmg_sw_taumol.f90` | Calculation of optical depths and Planck fractions for each spectral band
| `rrtmg_sw_vrtqdr.f90` | Two-stream vertical quadrature 
| `mcica_random_numbers.f90` | Random number generator for McICA
| `mcica_subcol_gen_sw.1col.f90` | Sub-column generator for McICA
| `rrtatm.f` | Process user-defined input data files
| `extra.f` | Process input data files
| `util_**.f` | Utilities (available for multiple platforms)

The following module files (in the `modules` directory) must be used to run RRTMG\_SW in stand-alone mode as a column model (these must be compiled before the source code files):

| File Name | Description |
| :--- | :--- |
| `parkind.f90` | real and integer kind type parameters |
| `parrrsw.f90` | main configuration parameters |
| `rrsw_aer.f90` | aerosol property coefficients |
| `rrsw_cld.f90` | cloud property coefficients |
| `rrsw_con.f90` | constants |
| `rrsw_kg**.f90` | absorption coefficient arrays for 16 spectral bands |
| `rrsw_ncpar.f90` | parameters for netCDF input data option |
| `rrsw_ref.f90` | reference atmosphere data arrays |
| `rrsw_tbl.f90` | exponential lookup table arrays |
| `rrsw_vsn.f90` | version number information |
| `rrsw_wvn.f90` | spectral band and *g*-interval array information |

### INPUT DATA
The following file (in the `data` directory) is the optional netCDF input file containing absorption coefficient and other input data for the model. The file is used if keyword `KGSRC` is set for netCDF input in the makefile. 

| File Name | Description |
| :--- | :--- |
| `rrtmg_sw.nc` | Optional netCDF input data file |

### MAKEFILES
The following files (in `build/makefiles` directory) can be used to compile RRTMG\_SW in stand-alone mode as a column model on various platforms.  Link one of these into the `build` directory to compile.

| File Name | Description |
| :--- | :--- |
| `make_rrtmg_sw_sgi` | Sample makefile for SGI
| `make_rrtmg_sw_sun` | Sample makefile for SUN
| `make_rrtmg_sw_linux_pgi` | Sample makefile for LINUX (PGI compiler)
| `make_rrtmg_sw_aix_xlf90` | Sample makefile for AIX (XLF90 compiler)
| `make_rrtmg_sw_OS_X_g95` | Sample makefile for OS_X (G95 compiler)
| `make_rrtmg_sw_OS_X_ibm_xl` | Sample makefile for OS_X (IBM XL compiler)

### SAMPLE INPUT/OUTPUT
Several sample input (and output) files are included in the `run_examples_std_atm` directory. Note that user-defined profiles may be used for as many as 200 layers.

| File Name | Description |
| :--- | :--- |
| `INPUT_RRTM` | Required input file for (clear sky) atmospheric specification |
| `IN_CLD_RRTM` | Required input file for cloud specification if clouds are present |
| `IN_AER_RRTM` |  Required input file for aerosol specification if aerosols are present |
| `OUTPUT_RRTM` | Main output file for atmospheric fluxes and heating rates |
| `input_rrtm.MLS-clr` | Sample 51 layer mid-latitude summer standard atmosphere |
| `input_rrtm.MLS-cld-imca0-icld2` | Sample 51 layer mid-latitude summer standard atmosphere with cloud flag turned on and maximum-random cloud overlap selected (without McICA) |
| `input_rrtm.MLS-cld-imca1-icld2` | Sample 51 layer mid-latitude summer standard atmosphere with cloud flag turned on and maximum-random cloud overlap selected (with McICA) |
| `input_rrtm.MLS-cld-imca1-icld5-idcor0` | Sample 51 layer mid-latitude summer standard atmosphere with cloud flag turned on and exponential-random cloud overlap and constant decorrelation length selected (with McICA) |
| `input_rrtm.MLS-cld-imca1-icld5-idcor1` | Sample 51 layer mid-latitude summer standard atmosphere with cloud flag turned on and exponential-random cloud overlap and varying decorrelation length selected (with McICA) |
| `input_rrtm.MLS-clr-aer12` | Sample 51 layer mid-latitude summer standard atmosphere with aerosol flag set |
| `input_rrtm.MLS-clr-sza45-isolvar0_tsi_avg` | Sample 51 layer mid-latitude summer standard atmosphere with solar zenith angle set to 45 degrees and using the NRLSSI2 solar source function with total solar irradiance for the mean solar cycle with no solar variability |
| `input_rrtm.MLS-clr-sza45-isolvar1_tsi_max` | Sample 51 layer mid-latitude summer standard atmosphere with solar zenith angle set to 45 degrees and using the NRLSSI2 solar source function with solar variability active and with total solar irradiance near the maximum in the mean solar cycle |
| `input_rrtm.MLS-clr-sza45-isolvar1_tsi_min` | Sample 51 layer mid-latitude summer standard atmosphere with solar zenith angle set to 45 degrees and using the  NRLSSI2 solar source function with solar variability active and with total solar irradiance near the minimum in the mean solar cycle |
| `input_rrtm.MLS-clr-sza45-isolvar2_01Jan1950` | Sample 51 layer mid-latitude summer standard atmosphere with solar zenith angle set to 45 degrees and using the NRLSSI2 solar source function with solar variability active and with total solar irradiance specified with facular and sunspot indices for 1 January 1950 |
| `input_rrtm.MLS-clr-sza45-isolvar3_bndscl_tsi_max` | Sample 51 layer mid-latitude summer standard atmosphere with solar zenith angle set to 45 degrees and using the NRLSSI2 solar source function with solar variability active and with total solar irradiance near the maximum in the mean solar cycle scaled to a different value with individual band scaling factors |
| `input_rrtm.MLW-clr` | Sample 51 layer mid-latitude winter standard atmosphere |
| `input_rrtm.SAW-clr` | Sample 51 layer sub-arctic winter standard atmosphere |
| `input_rrtm.TROP-clr` | Sample 51 layer tropical standard atmosphere |
| `in_cld_rrtm-cld5` | Sample cloud input file |
| `in_cld_rrtm-cld6` | Sample cloud input file |
| `in_cld_rrtm-cld7` | Sample cloud input file |
| `in_aer_rrtm-aer12` | Sample aerosol input file |
| `script.run_std_atm` | UNIX script for running the full suite of example cases, which will put the output into similarly named files prefixed with `output_rrtm*` |

### INSTRUCTIONS FOR COMPILING AND RUNNING THE COLUMN MODEL:
1) In the `build` directory, link one of the makefiles from the `makefile` sub-directory into `build/make.build`. To use the optional netCDF input file, switch the keyword `KGSRC` in the makefile from `dat` to `nc`. Compile the model with `make -f make.build`
2) Link the executable to, for example, `rrtmg_sw` in the `run_examples_std_atm` directory
3) If the optional netCDF input file was selected when compiling, link the file `data/rrtmg_sw.nc` into the `run_examples_std_atm` directory.  
4) In the `run_examples_std_atm` directory, run the UNIX script `./script.run_std_atm` to run the full suite of example cases. To run a single case, modify `INPUT_RRTM` following the instructions in `doc/rrtmg_sw_instructions.txt`, or copy one of the example `input_rrtm*` files into `INPUT_RRTM`. If clouds are selected (`ICLD` > 0), then modify `IN_CLD_RRTM` or copy one of the `in_cld_rrtm*` files into `IN_CLD_RRTM`. If aerosols are selected (`IAER` > 0), then modify `IN_AER_RRTM` or set it to the sample file `in_aer_rrtm-aer12`.
5) In column mode, if McICA is selected (`IMCA`=1) with partial cloudiness defined, then RRTMG\_SW will run the case 200 times to derive adequate statistics, and the average of the 200 samples will be written to the output file, `OUTPUT_RRTM`. 

## RRTMG_SW : GCM version <a name="gcm"></a>

### SOURCE CODE
The following source files (in the `src` directory) must be used to run RRTMG\_SW as a callable subroutine:

| File Name | Description |
| :--- | :--- |
| `rrtmg_sw_rad.f90` | RRTMG_SW main module (with McICA)
| `rrtmg_sw_rad.nomcica.f90` | Optional RRTMG_SW main module (without McICA only)
| `rrtmg_sw_cldprop.f90` | Calculation of cloud optical properties
| `rrtmg_sw_cldprmc.f90` | Calculation of cloud optical properties (McICA)
| `rrtmg_sw_init.f90` | RRTMG_SW initialization routine; reduces *g*-intervals from 224 to 112; (This has to run only once and should be installed in the GCM initialization section)
| `rrtmg_sw_k_g.f90` | Absorption coefficient data file
| `rrtmg_sw_read_nc.f90` | Alternate absorption coefficient data netCDF input
| `rrtmg_sw_reftra.f90` | Calculation of two-stream reflectivities and transmissivities
| `rrtmg_sw_setcoef.f90` | Set up routine
| `rrtmg_sw_spcvrt.f90` | Top subroutine for two-stream model
| `rrtmg_sw_spcvmc.f90` | Top subroutine for two-stream model (McICA)
| `rrtmg_sw_taumol.f90` | Calculation of optical depths and Planck fractions for each spectral band
| `rrtmg_sw_vrtqdr.f90` | Two-stream vertical quadrature 
| `mcica_random_numbers.f90` | Random number generator for McICA
| `mcica_subcol_gen_sw.f90` | Sub-column generator for McICA

**NOTE**: Only one of `rrtmg_sw_k_g.f90` or `rrtmg_sw_read_nc.f90` is required. 

The following module files (in the `modules` directory) must be used to run RRTMG\_SW as a callable subroutine (these must be compiled before the source code) 

| File Name | Description |
| :--- | :--- |
| `parkind.f90` | real and integer kind type parameters |
| `parrrsw.f90` | main configuration parameters |
| `rrsw_aer.f90` | aerosol property coefficients |
| `rrsw_cld.f90` | cloud property coefficients |
| `rrsw_con.f90` | constants |
| `rrsw_kg**.f90` | absorption coefficient arrays for 16 spectral bands |
| `rrsw_ncpar.f90` | parameters for netCDF input data option |
| `rrsw_ref.f90` | reference atmosphere data arrays |
| `rrsw_tbl.f90` | exponential lookup table arrays |
| `rrsw_vsn.f90` | version number information |
| `rrsw_wvn.f90` | spectral band and *g*-interval array information |

### INPUT DATA
The following file (in the `data` directory) is the optional netCDF file containing absorption coefficient and other input data for the model. The file is used if source file `rrtmg_sw_read_nc.f90` is used in place of `rrtmg_sw_k_g.f90` (only one or the other is required). 

| File Name | Description |
| :--- | :--- |
| `rrtmg_sw.nc` | Optional netCDF input data file |

### NOTES ON RUNNING THE GCM (SUBROUTINE) VERSION OF THE CODE

1) The module `rrtmg_sw_init.f90` is the initialization routine that has to be called only once.  The call to this subroutine should be moved to the initialization section of the host model if RRTMG_SW is called by a GCM or SCM. 
2) The number of model layers and the number of columns to be looped over should be passed into RRTMG_SW through the subroutine call along with the other model profile arrays.  
3) To utilize McICA, the sub-column generator (`mcica_subcol_gen_sw.f90`) must be implemented in the GCM so that it is called just before RRTMG\_SW. The cloud overlap method is selected using the input flag, icld. If either exponential (`ICLD`=4) or exponential-random (`ICLD`=5) cloud overlap is selected, then the  subroutine `get_alpha` must be called prior to calling `mcica_subcol_sw` to define the vertical correlation parameter, `alpha`, needed for those overlap methods. Also for those methods, use the input flag `idcor` to select the use of either a constant or latitude-varying decorrelation length. If McICA is utilized, this will run only a single statistical sample per model grid box. There are two options for the random number generator used with McICA, which is selected with the variable `irnd` in `mcica_subcol_gen_sw.f90`. When using McICA, then the main module is `rrtmg_sw_rad.f90`. If McICA is not used, then the main module is `rrtmg_sw_rad.nomcica.f90`, though the cloud specification is limited to overcast clouds.

## Maintenance and Contact Info <a name="contact"></a>

Atmospheric and Environmental Research 
131 Hartwell Avenue, Lexington, MA 02421

Original version:   Eli. J. Mlawer, J. S. Delamere, et al. (AER)
Revision for GCMs:  Michael J. Iacono (AER)

Contact:   Michael J. Iacono   (E-mail: miacono@aer.com)

## References <a name="ref"></a>

* [AER Radiative Transfer Models Documentation](https://www.rtweb.aer.com)
* [Github Wiki](https://github.com/AER-RC/RRTMG_LW/wiki)
* **RRTMG_SW, RRTM_SW**
  * Iacono, M.J., J.S. Delamere, E.J. Mlawer, M.W. Shephard, S.A. Clough, and W.D. Collins, Radiative forcing by long-lived greenhouse gases: Calculations with the AER radiative transfer models, *J. Geophys. Res.*, 113, D13103, doi:10.1029/2008JD009944, 2008.
  * Clough, S.A., M.W. Shephard, E.J. Mlawer, J.S. Delamere, M.J. Iacono, K. Cady-Pereira, S. Boukabara, and P.D. Brown, Atmospheric radiative transfer modeling: a summary of the AER codes, *J. Quant. Spectrosc. Radiat. Transfer*, 91, 233-244, 2005. 
* **McICA**
  * Pincus, R., H. W. Barker, and J.-J. Morcrette, A fast, flexible, approximation technique for computing radiative transfer in inhomogeneous cloud fields, *J. Geophys. Res.*, 108(D13), 4376, doi:10.1029/2002JD003322, 2003.
*  **Latitude-Varying Decorrelation Length**
    *  Oreopoulos, L., D. Lee, Y.C. Sud, and M.J. Suarez, Radiative impacts of cloud heterogeneity and overlap in an atmospheric General Circulation Model, *Atmos. Chem. Phys.*, 12, 9097-9111, doi:10.5194/acp-12-9097-2012, 2012.
* [Full list of references](https://github.com/AER-RC/RRTMG_SW/wiki/References)
