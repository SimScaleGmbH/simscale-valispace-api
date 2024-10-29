from typing import Any, Dict

# Simscale Libraries
from simscale_sdk import Configuration, ApiClient, ProjectsApi, StorageApi, GeometryImportsApi, GeometriesApi, MeshOperationsApi, SimulationsApi, SimulationRunsApi, Project, GeometryImportRequest, ApiException, ReportsApi, MaterialsApi
from simscale_sdk import CoupledConjugateHeatTransfer, FluidModel, DimensionalVectorAcceleration, FluidInitialConditions, AdvancedConcepts, ConvectiveHeatTransferMaterials, TopologicalReference, FluidNumerics, RelaxationFactor, DimensionalPressure, ResidualControls, Tolerance, FluidSolvers, Schemes, TimeDifferentiationSchemes, GradientSchemes, DivergenceSchemes, LaplacianSchemes, InterpolationSchemes, SurfaceNormalGradientSchemes, VelocityInletBC, FixedValueVBC, DimensionalVectorFunctionSpeed, ComponentVectorFunction, ConstantFunction, FixedValueTBC, DimensionalFunctionTemperature, PressureOutletBC, FixedValuePBC, DimensionalFunctionPressure, WallBC, NoSlipVBC, FluidSimulationControl, DimensionalTime, TimeStepWriteControl, ScotchDecomposeAlgorithm, FluidResultControls, AreaAverageResultControl, ProbePointsResultControl
from simscale_sdk import GeometryImportRequestLocation, GeometryImportRequestOptions, Point, DimensionalVectorLength, DecimalVector
from simscale_sdk import SimulationSpec, MeshOperation, SimmetrixMeshingFluid, AutomaticLayerOn, SimulationRun
from simscale_sdk import StationaryTimeDependency, FluidInterface, RegionInterface, CoupledInterfaceThermal,CoupledConjugateHeatTransferMaterials, IncompressibleMaterial, NewtonianViscosityModel,  IncompressibleFluidMaterials
from simscale_sdk import DimensionalKinematicViscosity, DimensionalDensity, DimensionalKinematicViscosity, DimensionalThermalExpansionRate
from simscale_sdk import DimensionalTemperature, DimensionalSpecificHeat, TopologicalReference, SolidCompressibleMaterial, DimensionalFunctionThermalConductivity, HConstThermo, DimensionalSpecificHeat, RhoConstEquationOfState
from simscale_sdk import ConstIsoTransport, IsotropicConductivity, DimensionalInitialConditionDomainsPressure, DimensionalPressure, DimensionalVectorInitialConditionDomainsSpeed, DimensionalVectorSpeed, DimensionalInitialConditionDomainsTemperature, DimensionalInitialConditionDomainsTurbulenceKineticEnergy, DimensionalTurbulenceKineticEnergy, DimensionalInitialConditionDomainsSpecificTurbulenceDissipationRate, DimensionalSpecificTurbulenceDissipationRate, FlowRateInletVBC, MassFlow, DimensionalFunctionMassFlowRate, ExternalWallHeatFluxTBC, FixedPowerHeatFlux, DimensionalFunctionPower, PBICGSolver, DILUPreconditioner, PBICGSolver, ILUCpPreconditioner, GAMGSolver, PBICGSolver, DILUPreconditioner
from simscale_sdk import MaterialLibraryReference, Stabilization, FieldLimits, AutomaticTurbulence, ConstAnIsoTransport
from simscale_sdk import OrthotropicConductivity, CartesianOrientation, DerivedHeatFlux, DimensionalThermalTransmittance, NoWallThermal
from simscale_sdk import AbsolutePowerSource, SmoothSolver, MeshesApi

import time
import urllib3
import csv
import valispace
import re



# Importing the SIMSCALE_API_KEY from user secrets defined in Settings
from .settings import simscale_key

# Importing Valispace username and password to be used in the script
from .settings import Username, Password

def main(**kwargs) -> Dict[str, Any]:
    """
    This is the main function to execute your script and it must exists.

    Other functions and files can be also created. You have at your disposal
    to import Valispace API, scipy, numpy and pint.

    To import user secrets use : from .settings import name_of_secret
        - Example from .settings import USERNAME, PASSWORD (the secrets need to be set in the user settings
    )
    :param kwargs: Dictionary with data received from Valispace.
    :type kwargs: Dict[str, Any]

    :return: Dictionary with data to send back to Valispace.
    :rtype: Dict[str, Any]
    """
    # Feed the SimScale API relevant input
    SIMSCALE_API_URL = "https://api.simscale.com"
    SIMSCALE_API_KEY = simscale_key

    # Confirm environment variables are assigned
    if not SIMSCALE_API_KEY or not SIMSCALE_API_URL:
        raise Exception("Either `SIMSCALE_API_KEY` or `SIMSCALE_API_URL` environment variable is missing. Add the API key to the text document key.txt")

    # API client configuration
    api_key_header = "X-API-KEY"
    api_key = SIMSCALE_API_KEY
    configuration = Configuration()
    configuration.host = SIMSCALE_API_URL + "/v0"
    configuration.api_key = {
        api_key_header: api_key,
    }
    configuration.debug = True

    api_client = ApiClient(configuration)

    retry_policy = urllib3.Retry(connect=5, read=5, redirect=0, status=5, backoff_factor=0.2)
    api_client.rest_client.pool_manager.connection_pool_kw["retries"] = retry_policy

    # SimScale API clients
    project_api = ProjectsApi(api_client)
    storage_api = StorageApi(api_client)
    geometry_import_api = GeometryImportsApi(api_client)
    geometry_api = GeometriesApi(api_client)
    mesh_operation_api = MeshOperationsApi(api_client)
    meshes_api = MeshesApi(api_client)
    simulation_api = SimulationsApi(api_client)
    simulation_run_api = SimulationRunsApi(api_client)
    #reports_api = ReportsApi(api_client)
    #materials_api = MaterialsApi(api_client)


    valispace_api = valispace.API(url="https://simscale.valispace.com/",session_token=kwargs['temporary_access_token'])

    #****** VALISPACE GET VALUES **********

    #Get Power, Flow rate and Material values from valispace
    cellPower_vali_ID = 12921
    cellPower  = valispace_api.get(f"valis/{cellPower_vali_ID}")
    cellPower_value = cellPower['value']

    massFlowRate_vali_ID = 12922
    massFlowRate = valispace_api.get(f"valis/{massFlowRate_vali_ID}") 
    massFlowRate_value = massFlowRate['value']

    materialName_vali_ID = 105
    materialCp_vali_ID = 12945
    materialDensity_vali_ID = 12946
    materialViscosity_vali_ID = 12947

    [materialName, materialCP, materialDensity, materialViscosity] = [valispace_api.get(f"textvalis/{materialName_vali_ID}"), 
                                                    valispace_api.get(f"valis/{materialCp_vali_ID}"), 
                                                    valispace_api.get(f"valis/{materialDensity_vali_ID}"), 
                                                    valispace_api.get(f"valis/{materialViscosity_vali_ID}") ]
    [materialName_value, materialCp_value, materialDensity_value, materialViscosity_value] = [materialName['text'], materialCP['value'], materialDensity['value'], materialViscosity['value']]
    #****** VALISPACE GET VALUES END **********

    
    #Convert to Simscale compatible values to be fed into the simulation specifications

    #**Drive Motor** 
    power_value = cellPower_value
    power = str(cellPower_value)+"W"
  

    #**Flow Profile**
    flowrate_value = massFlowRate_value
    flowrate = str(massFlowRate_value)+"kg/s"


    #**Material Definitation**
    pattern = r'<p>(.*?)<\/p>'
    converted_text = re.findall(pattern, materialName_value, re.DOTALL) #convert from html to plain text     
    material_name = converted_text[0]
    kinematic_viscosity = materialViscosity_value
    specific_heat = materialCp_value
    density = materialDensity_value


    #******** Set SimScale project, geometry, meshing ID. *************************
    
    # Hard-code these IDs for this project - This can be also looked up by its name. 
    # The project ID can be determined from the URL - check after "pid=" 
    project_id = "7984451083890346384" #Webinar project
    print(f"project_id: {project_id}")

    # Hard-coded geometry ID. - This can be also looked up by its name. 
    # Hint for now: Geometry ID can be retrieved through the browser console. This might change in the future.
    geometry_id = "fdeb71d3-8104-4330-87e9-52e4fbd7bb08"


    # Hard-coded meshing ID. 
    # Hint for now: The meshing ID can be retrieved through its name
    mesh_name = "Default Mesh"

    #******** Set SimScale project, geometry, meshing ID - END *************************



    # Read project information and update with the deserialized model
    project = project_api.get_project(project_id)
    project_api.update_project(project_id, project)

    # Geometry
    # TODO: Read existing geometries and pick the geometry_id from a given name
    # instead of hard-coding it 
    geometry_id = geometry_id
    print(f"geometry_id: {geometry_id}")

    # Read geometry information and update with the deserialized model (not used atm)
    geometry = geometry_api.get_geometry(project_id, geometry_id)
    geometry_api.update_geometry(project_id, geometry_id, geometry)


    # Create simulation specs - CHTv2 (Conjugate Heat Transfer)
    # IMPORTANT: This part was created using the get_sdk_code_example.py script except
    #            from the point monitoring part which was altered a bit.
    #
    # This part cotnains all the relevent simulation information e.g. Boundary conditions, numerics, Result controls etc

    model = CoupledConjugateHeatTransfer(
        connection_groups=[],
        model=FluidModel(
            gravity=DimensionalVectorAcceleration(
                value=DecimalVector(
                    x=0,
                    y=0,
                    z=0,
                ),
                unit="m/s²",
            ),
        ),
        materials=CoupledConjugateHeatTransferMaterials(
            fluids=[
                IncompressibleMaterial(
                    type="INCOMPRESSIBLE",
                    name=material_name,
                    viscosity_model=NewtonianViscosityModel(
                        type="NEWTONIAN",
                        kinematic_viscosity=DimensionalKinematicViscosity(
                            value=kinematic_viscosity,
                            unit="m²/s",
                        ),
                    ),
                    density=DimensionalDensity(
                        value=density,
                        unit="kg/m³",
                    ),
                    thermal_expansion_coefficient=DimensionalThermalExpansionRate(
                        value=2.07E-4,
                        unit="1/K",
                    ),
                    reference_temperature=DimensionalTemperature(
                        value=298.15,
                        unit="K",
                    ),
                    laminar_prandtl_number=6.5241,
                    turbulent_prandtl_number=0.85,
                    specific_heat=DimensionalSpecificHeat(
                        value=specific_heat,
                        unit="J/(kg·K)",
                    ),
                    topological_reference=TopologicalReference(
                        entities=[
                            "B11_TE1",
                            "B3_TE1",
                        ],
                        sets=[],
                    ),
                    material_library_reference=MaterialLibraryReference(
                        material_group_id="internal:af500b1b-bafb-421d-afa3-c976dc706f52",
                        material_id="internal:f6c40e2d-b832-473a-b0de-b8b732e50a10",
                        interpolation_parameters={},
                    ),
                ),
            ],
            solids=[
                SolidCompressibleMaterial(
                    name="Aluminium",
                    transport=ConstIsoTransport(
                        type="CONST_ISO",
                        conductivity=IsotropicConductivity(
                            type="ISOTROPIC",
                            thermal_conductivity=DimensionalFunctionThermalConductivity(
                                value=ConstantFunction(
                                    type="CONSTANT",
                                    value=235,
                                ),
                                unit="W/(m·K)",
                            ),
                        ),
                        thermo=HConstThermo(
                            type="HCONST",
                            specific_heat=DimensionalSpecificHeat(
                                value=897,
                                unit="J/(kg·K)",
                            ),
                            equation_of_state=RhoConstEquationOfState(
                                type="RHO_CONST",
                                density=DimensionalDensity(
                                    value=2700,
                                    unit="kg/m³",
                                ),
                            ),
                        ),
                    ),
                    topological_reference=TopologicalReference(
                        entities=[
                            "B16_TE138",
                            "B4_TE138",
                        ],
                        sets=[],
                    ),
                    material_library_reference=MaterialLibraryReference(
                        material_group_id="internal:af500b1b-bafb-421d-afa3-c976dc706f52",
                        material_id="internal:0cd6d090-6cb3-4488-a3b1-ef7838f38923",
                        interpolation_parameters={},
                    ),
                ),
                SolidCompressibleMaterial(
                    name="Copper",
                    transport=ConstIsoTransport(
                        type="CONST_ISO",
                        conductivity=IsotropicConductivity(
                            type="ISOTROPIC",
                            thermal_conductivity=DimensionalFunctionThermalConductivity(
                                value=ConstantFunction(
                                    type="CONSTANT",
                                    value=401,
                                ),
                                unit="W/(m·K)",
                            ),
                        ),
                        thermo=HConstThermo(
                            type="HCONST",
                            specific_heat=DimensionalSpecificHeat(
                                value=385,
                                unit="J/(kg·K)",
                            ),
                            equation_of_state=RhoConstEquationOfState(
                                type="RHO_CONST",
                                density=DimensionalDensity(
                                    value=8960,
                                    unit="kg/m³",
                                ),
                            ),
                        ),
                    ),
                    topological_reference=TopologicalReference(
                        entities=[
                            "B104_TE36251",
                            "B10_TE36251",
                            "B7_TE22976",
                            "B8_TE56196",
                            "B6_TE82777",
                            "B80_TE62583",
                            "B121_TE15464",
                            "B5_TE62583",
                            "B115_TE82777",
                            "B157_TE56196",
                            "B38_TE22976",
                            "B9_TE15464",
                        ],
                        sets=[],
                    ),
                    material_library_reference=MaterialLibraryReference(
                        material_group_id="internal:af500b1b-bafb-421d-afa3-c976dc706f52",
                        material_id="internal:e9c82c2b-339a-4bf6-a5b8-a9dcb74b9e58",
                        interpolation_parameters={},
                    ),
                ),
                SolidCompressibleMaterial(
                    name="Battery Cells",
                    transport=ConstAnIsoTransport(
                        type="CONST_AN_ISO",
                        conductivity=OrthotropicConductivity(
                            type="ORTHOTROPIC",
                            thermal_conductivity_x=DimensionalFunctionThermalConductivity(
                                value=ConstantFunction(
                                    type="CONSTANT",
                                    value=0.6,
                                ),
                                unit="W/(m·K)",
                            ),
                            thermal_conductivity_y=DimensionalFunctionThermalConductivity(
                                value=ConstantFunction(
                                    type="CONSTANT",
                                    value=0.6,
                                ),
                                unit="W/(m·K)",
                            ),
                            thermal_conductivity_z=DimensionalFunctionThermalConductivity(
                                value=ConstantFunction(
                                    type="CONSTANT",
                                    value=40,
                                ),
                                unit="W/(m·K)",
                            ),
                        ),
                        orientation=CartesianOrientation(
                            type="CARTESIAN",
                        ),
                        thermo=HConstThermo(
                            type="HCONST",
                            specific_heat=DimensionalSpecificHeat(
                                value=1030,
                                unit="J/(kg·K)",
                            ),
                            equation_of_state=RhoConstEquationOfState(
                                type="RHO_CONST",
                                density=DimensionalDensity(
                                    value=2300,
                                    unit="kg/m³",
                                ),
                            ),
                        ),
                    ),
                    topological_reference=TopologicalReference(
                        entities=[
                            "B61_TE53580",
                            "B169_TE13400",
                            "B171_TE32829",
                            "B151_TE52516",
                            "B161_TE12124",
                            "B42_TE12703",
                            "B144_TE60330",
                            "B143_TE79414",
                            "B81_TE32216",
                            "B178_TE78947",
                            "B122_TE78879",
                            "B58_TE79082",
                            "B139_TE80029",
                            "B74_TE59461",
                            "B13_TE32891",
                            "B170_TE19461",
                            "B138_TE19817",
                            "B136_TE33017",
                            "B88_TE18661",
                            "B87_TE12065",
                            "B30_TE11948",
                            "B37_TE59768",
                            "B46_TE60112",
                            "B35_TE59811",
                            "B130_TE32173",
                            "B59_TE79601",
                            "B114_TE19588",
                            "B127_TE12901",
                            "B159_TE53535",
                            "B55_TE32000",
                            "B118_TE59854",
                            "B117_TE31871",
                            "B89_TE59229",
                            "B154_TE18929",
                            "B124_TE79855",
                            "B134_TE18448",
                            "B146_TE60026",
                            "B17_TE11370",
                            "B77_TE12301",
                            "B28_TE59639",
                            "B168_TE32346",
                            "B65_TE78748",
                            "B48_TE32955",
                            "B54_TE19274",
                            "B56_TE52985",
                            "B62_TE53737",
                            "B78_TE32086",
                            "B63_TE53274",
                            "B109_TE11659",
                            "B140_TE52929",
                            "B44_TE18861",
                            "B123_TE59553",
                            "B99_TE19648",
                            "B53_TE32648",
                            "B57_TE53682",
                            "B15_TE52732",
                            "B72_TE18590",
                            "B153_TE19931",
                            "B165_TE79351",
                            "B173_TE32043",
                            "B68_TE52646",
                            "B49_TE13340",
                            "B76_TE11835",
                            "B126_TE53362",
                            "B128_TE59897",
                            "B106_TE13444",
                            "B162_TE18724",
                            "B107_TE79016",
                            "B51_TE31957",
                            "B33_TE59725",
                            "B82_TE60373",
                            "B108_TE19400",
                            "B174_TE31914",
                            "B26_TE79539",
                            "B90_TE52689",
                            "B92_TE53491",
                            "B131_TE32734",
                            "B148_TE11717",
                            "B142_TE12008",
                            "B85_TE12491",
                            "B50_TE59682",
                            "B83_TE52472",
                            "B69_TE11317",
                            "B163_TE79226",
                            "B73_TE78475",
                            "B84_TE60244",
                            "B116_TE18376",
                            "B12_TE19993",
                            "B160_TE32777",
                            "B47_TE59369",
                            "B181_TE11891",
                            "B94_TE18515",
                            "B36_TE79665",
                            "B97_TE59323",
                            "B67_TE11599",
                            "B105_TE53405",
                            "B43_TE52874",
                            "B141_TE11427",
                            "B24_TE32302",
                            "B93_TE59276",
                            "B156_TE32432",
                            "B113_TE79973",
                            "B102_TE60200",
                            "B149_TE19005",
                            "B19_TE78543",
                            "B166_TE32691",
                            "B79_TE32129",
                            "B176_TE53318",
                            "B100_TE59415",
                            "B172_TE19525",
                            "B111_TE79791",
                            "B32_TE60155",
                            "B39_TE19074",
                            "B175_TE53448",
                            "B150_TE19705",
                            "B22_TE53230",
                            "B98_TE13248",
                            "B177_TE53177",
                            "B125_TE79916",
                            "B41_TE18794",
                            "B31_TE18305",
                            "B155_TE78409",
                            "B101_TE59596",
                            "B23_TE52818",
                            "B96_TE52775",
                            "B18_TE19140",
                            "B167_TE59509",
                            "B27_TE59940",
                            "B132_TE53052",
                            "B60_TE79154",
                            "B152_TE11778",
                            "B71_TE11542",
                            "B52_TE19760",
                            "B70_TE32389",
                            "B179_TE32562",
                            "B20_TE52560",
                            "B34_TE52603",
                            "B66_TE13160",
                            "B40_TE79476",
                            "B119_TE80086",
                            "B112_TE79289",
                            "B110_TE60069",
                            "B147_TE32519",
                            "B29_TE13068",
                            "B21_TE32259",
                            "B137_TE12242",
                            "B158_TE32475",
                            "B164_TE60287",
                            "B86_TE53110",
                            "B45_TE12182",
                            "B133_TE53623",
                            "B75_TE32605",
                            "B64_TE78813",
                            "B120_TE19205",
                            "B145_TE33081",
                            "B95_TE79727",
                            "B91_TE59983",
                            "B135_TE78606",
                            "B103_TE11485",
                            "B25_TE19873",
                            "B129_TE19338",
                            "B180_TE78671",
                        ],
                        sets=[],
                    ),
                    material_library_reference=MaterialLibraryReference(
                        material_group_id="internal:af500b1b-bafb-421d-afa3-c976dc706f52",
                        material_id="internal:baf621d5-4ee5-4fe6-beb5-a36813cdea40",
                        interpolation_parameters={},
                    ),
                ),
                SolidCompressibleMaterial(
                    name="Polypropylene",
                    transport=ConstIsoTransport(
                        type="CONST_ISO",
                        conductivity=IsotropicConductivity(
                            type="ISOTROPIC",
                            thermal_conductivity=DimensionalFunctionThermalConductivity(
                                value=ConstantFunction(
                                    type="CONSTANT",
                                    value=0.2,
                                ),
                                unit="W/(m·K)",
                            ),
                        ),
                        thermo=HConstThermo(
                            type="HCONST",
                            specific_heat=DimensionalSpecificHeat(
                                value=1920,
                                unit="J/(kg·K)",
                            ),
                            equation_of_state=RhoConstEquationOfState(
                                type="RHO_CONST",
                                density=DimensionalDensity(
                                    value=900,
                                    unit="kg/m³",
                                ),
                            ),
                        ),
                    ),
                    topological_reference=TopologicalReference(
                        entities=[
                            "B1_TE60416",
                            "B14_TE60416",
                        ],
                        sets=[],
                    ),
                    material_library_reference=MaterialLibraryReference(
                        material_group_id="internal:af500b1b-bafb-421d-afa3-c976dc706f52",
                        material_id="internal:951db993-29d8-41b8-874b-4fc209daf402",
                        interpolation_parameters={},
                    ),
                ),
                SolidCompressibleMaterial(
                    name="Potting",
                    transport=ConstIsoTransport(
                        type="CONST_ISO",
                        conductivity=IsotropicConductivity(
                            type="ISOTROPIC",
                            thermal_conductivity=DimensionalFunctionThermalConductivity(
                                value=ConstantFunction(
                                    type="CONSTANT",
                                    value=20,
                                ),
                                unit="W/(m·K)",
                            ),
                        ),
                        thermo=HConstThermo(
                            type="HCONST",
                            specific_heat=DimensionalSpecificHeat(
                                value=450,
                                unit="J/(kg·K)",
                            ),
                            equation_of_state=RhoConstEquationOfState(
                                type="RHO_CONST",
                                density=DimensionalDensity(
                                    value=795,
                                    unit="kg/m³",
                                ),
                            ),
                        ),
                    ),
                    topological_reference=TopologicalReference(
                        entities=[
                            "B2_TE125",
                        ],
                        sets=[],
                    ),
                    material_library_reference=MaterialLibraryReference(
                        material_group_id="internal:af500b1b-bafb-421d-afa3-c976dc706f52",
                        material_id="internal:951db993-29d8-41b8-874b-4fc209daf402",
                        interpolation_parameters={},
                    ),
                ),
            ],
        ),
        initial_conditions=FluidInitialConditions(
            gauge_pressure_rgh=DimensionalInitialConditionDomainsPressure(
                _global=DimensionalPressure(
                    value=0,
                    unit="Pa",
                ),
                subdomains=[],
            ),
            velocity=DimensionalVectorInitialConditionDomainsSpeed(
                _global=DimensionalVectorSpeed(
                    value=DecimalVector(
                        x=0,
                        y=0,
                        z=0,
                    ),
                    unit="m/s",
                ),
                subdomains=[],
            ),
            temperature=DimensionalInitialConditionDomainsTemperature(
                _global=DimensionalTemperature(
                    value=19.85,
                    unit="°C",
                ),
                subdomains=[],
            ),
            turbulent_kinetic_energy=DimensionalInitialConditionDomainsTurbulenceKineticEnergy(
                _global=DimensionalTurbulenceKineticEnergy(
                    value=0.00375,
                    unit="m²/s²",
                ),
            ),
            omega_dissipation_rate=DimensionalInitialConditionDomainsSpecificTurbulenceDissipationRate(
                _global=DimensionalSpecificTurbulenceDissipationRate(
                    value=3.375,
                    unit="1/s",
                ),
            ),
        ),
        boundary_conditions=[
            VelocityInletBC(
                name="Velocity inlet 1",
                type="VELOCITY_INLET_V3",
                velocity=FlowRateInletVBC(
                    type="FLOW_RATE_INLET_VELOCITY",
                    flow_rate=MassFlow(
                        type="MASS",
                        value=DimensionalFunctionMassFlowRate(
                            value=ConstantFunction(
                                type="CONSTANT",
                                value=flowrate_value,
                            ),
                            unit="kg/s",
                        ),
                    ),
                ),
                turbulence=AutomaticTurbulence(
                    type="AUTOMATIC_TURBULENCE",
                ),
                temperature=FixedValueTBC(
                    type="FIXED_VALUE",
                    value=DimensionalFunctionTemperature(
                        value=ConstantFunction(
                            type="CONSTANT",
                            value=18,
                        ),
                        unit="°C",
                    ),
                ),
                topological_reference=TopologicalReference(
                    entities=[
                        "B3_TE1883",
                        "B11_TE1883",
                    ],
                    sets=[],
                ),
            ),
            PressureOutletBC(
                name="Pressure outlet 2",
                type="PRESSURE_OUTLET_V30",
                gauge_pressure_rgh=FixedValuePBC(
                    type="FIXED_VALUE",
                    value=DimensionalFunctionPressure(
                        value=ConstantFunction(
                            type="CONSTANT",
                            value=0,
                        ),
                        unit="Pa",
                    ),
                ),
                topological_reference=TopologicalReference(
                    entities=[
                        "B3_TE4396",
                        "B11_TE4396",
                    ],
                    sets=[],
                ),
            ),
            WallBC(
                name="External Walls",
                type="WALL_V34",
                velocity=NoSlipVBC(
                    type="NO_SLIP",
                ),
                temperature=ExternalWallHeatFluxTBC(
                    type="EXTERNAL_WALL_HEAT_FLUX_TEMPERATURE",
                    heat_flux=DerivedHeatFlux(
                        type="DERIVED",
                        heat_transfer_coefficient=DimensionalThermalTransmittance(
                            value=5,
                            unit="W/(K·m²)",
                        ),
                        ambient_temperature=DimensionalTemperature(
                            value=30,
                            unit="°C",
                        ),
                        wall_thermal=NoWallThermal(
                            type="NO_RESISTANCE",
                        ),
                    ),
                ),
                topological_reference=TopologicalReference(
                    entities=[
                        "B4_TE687",
                        "B16_TE687",
                        "B16_TE304",
                        "B16_TE300",
                        "B16_TE10648",
                        "B16_TE68",
                        "B16_TE308",
                        "B16_TE426",
                        "B16_TE448",
                        "B16_TE10578",
                        "B16_TE12684",
                        "B16_TE122",
                        "B16_TE12674",
                        "B16_TE12679",
                        "B16_TE1",
                        "B4_TE426",
                        "B4_TE122",
                        "B4_TE12674",
                        "B4_TE12679",
                        "B4_TE12684",
                        "B4_TE300",
                        "B4_TE10648",
                        "B4_TE308",
                        "B4_TE10578",
                        "B4_TE68",
                        "B4_TE1",
                        "B4_TE448",
                        "B4_TE304",
                        "B2_TE212",
                        "B2_TE121",
                        "B2_TE206",
                        "B2_TE209",
                    ],
                    sets=[],
                ),
            ),
        ],
        advanced_concepts=AdvancedConcepts(
            porous_mediums=[],
            power_sources=[
                AbsolutePowerSource(
                    type="ABSOLUTE_V23",
                    name="Absolute power source 1",
                    heat_flux=DimensionalFunctionPower(
                        value=ConstantFunction(
                            type="CONSTANT",
                            value=power_value,
                        ),
                        unit="W",
                    ),
                    topological_reference=TopologicalReference(
                        entities=[
                            "B61_TE53580",
                            "B169_TE13400",
                            "B171_TE32829",
                            "B151_TE52516",
                            "B161_TE12124",
                            "B42_TE12703",
                            "B144_TE60330",
                            "B143_TE79414",
                            "B81_TE32216",
                            "B178_TE78947",
                            "B122_TE78879",
                            "B58_TE79082",
                            "B139_TE80029",
                            "B74_TE59461",
                            "B13_TE32891",
                            "B170_TE19461",
                            "B138_TE19817",
                            "B136_TE33017",
                            "B88_TE18661",
                            "B87_TE12065",
                            "B30_TE11948",
                            "B37_TE59768",
                            "B46_TE60112",
                            "B35_TE59811",
                            "B130_TE32173",
                            "B59_TE79601",
                            "B114_TE19588",
                            "B127_TE12901",
                            "B159_TE53535",
                            "B55_TE32000",
                            "B118_TE59854",
                            "B117_TE31871",
                            "B89_TE59229",
                            "B154_TE18929",
                            "B124_TE79855",
                            "B134_TE18448",
                            "B146_TE60026",
                            "B17_TE11370",
                            "B77_TE12301",
                            "B28_TE59639",
                            "B168_TE32346",
                            "B65_TE78748",
                            "B48_TE32955",
                            "B54_TE19274",
                            "B56_TE52985",
                            "B62_TE53737",
                            "B78_TE32086",
                            "B63_TE53274",
                            "B109_TE11659",
                            "B140_TE52929",
                            "B44_TE18861",
                            "B123_TE59553",
                            "B99_TE19648",
                            "B53_TE32648",
                            "B57_TE53682",
                            "B15_TE52732",
                            "B72_TE18590",
                            "B153_TE19931",
                            "B165_TE79351",
                            "B173_TE32043",
                            "B68_TE52646",
                            "B49_TE13340",
                            "B76_TE11835",
                            "B126_TE53362",
                            "B128_TE59897",
                            "B106_TE13444",
                            "B162_TE18724",
                            "B107_TE79016",
                            "B51_TE31957",
                            "B33_TE59725",
                            "B82_TE60373",
                            "B108_TE19400",
                            "B174_TE31914",
                            "B26_TE79539",
                            "B90_TE52689",
                            "B92_TE53491",
                            "B131_TE32734",
                            "B148_TE11717",
                            "B142_TE12008",
                            "B85_TE12491",
                            "B50_TE59682",
                            "B83_TE52472",
                            "B69_TE11317",
                            "B163_TE79226",
                            "B73_TE78475",
                            "B84_TE60244",
                            "B116_TE18376",
                            "B12_TE19993",
                            "B160_TE32777",
                            "B47_TE59369",
                            "B181_TE11891",
                            "B94_TE18515",
                            "B36_TE79665",
                            "B97_TE59323",
                            "B67_TE11599",
                            "B105_TE53405",
                            "B43_TE52874",
                            "B141_TE11427",
                            "B24_TE32302",
                            "B93_TE59276",
                            "B156_TE32432",
                            "B113_TE79973",
                            "B102_TE60200",
                            "B149_TE19005",
                            "B19_TE78543",
                            "B166_TE32691",
                            "B79_TE32129",
                            "B176_TE53318",
                            "B100_TE59415",
                            "B172_TE19525",
                            "B111_TE79791",
                            "B32_TE60155",
                            "B39_TE19074",
                            "B175_TE53448",
                            "B150_TE19705",
                            "B22_TE53230",
                            "B98_TE13248",
                            "B177_TE53177",
                            "B125_TE79916",
                            "B41_TE18794",
                            "B31_TE18305",
                            "B155_TE78409",
                            "B101_TE59596",
                            "B23_TE52818",
                            "B96_TE52775",
                            "B18_TE19140",
                            "B167_TE59509",
                            "B27_TE59940",
                            "B132_TE53052",
                            "B60_TE79154",
                            "B152_TE11778",
                            "B71_TE11542",
                            "B52_TE19760",
                            "B70_TE32389",
                            "B179_TE32562",
                            "B20_TE52560",
                            "B34_TE52603",
                            "B66_TE13160",
                            "B40_TE79476",
                            "B119_TE80086",
                            "B112_TE79289",
                            "B110_TE60069",
                            "B147_TE32519",
                            "B29_TE13068",
                            "B21_TE32259",
                            "B137_TE12242",
                            "B158_TE32475",
                            "B164_TE60287",
                            "B86_TE53110",
                            "B45_TE12182",
                            "B133_TE53623",
                            "B75_TE32605",
                            "B64_TE78813",
                            "B120_TE19205",
                            "B145_TE33081",
                            "B95_TE79727",
                            "B91_TE59983",
                            "B135_TE78606",
                            "B103_TE11485",
                            "B25_TE19873",
                            "B129_TE19338",
                            "B180_TE78671",
                        ],
                        sets=[],
                    ),
                    geometry_primitive_uuids=[],
                ),
            ],
            momentum_sources=[],
            thermal_resistance_networks=[],
        ),
        numerics=FluidNumerics(
            relaxation_type="MANUAL",
            relaxation_factor=RelaxationFactor(
                pressure_rgh_field=0.7,
                velocity_equation=0.3,
                temperature_equation=0.8,
                turbulent_kinetic_energy_equation=0.3,
                omega_dissipation_rate_equation=0.3,
            ),
            num_non_orthogonal_correctors=1,
            solvers=FluidSolvers(
                velocity_solver=PBICGSolver(
                    type="PBICG",
                    absolute_tolerance=1.0E-8,
                    relative_tolerance=0.01,
                    preconditioner=DILUPreconditioner(
                        type="DILU",
                    ),
                ),
                temperature_solver=PBICGSolver(
                    type="PBICG",
                    absolute_tolerance=1.0E-6,
                    relative_tolerance=0.1,
                    preconditioner=ILUCpPreconditioner(
                        type="ILUCP",
                        fill_in_level=1,
                    ),
                ),
                pressure_rgh_solver=GAMGSolver(
                    type="GAMG",
                    absolute_tolerance=1.0E-8,
                    relative_tolerance=0.01,
                    smoother="GAUSSSEIDEL",
                    num_pre_sweeps=1,
                    num_post_sweeps=1,
                    cache_agglomeration_on=True,
                    num_cells_coarsest_level=100,
                    num_merge_levels=1,
                ),
                turbulent_kinetic_energy_solver=SmoothSolver(
                    type="SMOOTH",
                    absolute_tolerance=1.0E-8,
                    relative_tolerance=0.01,
                    smoother="GAUSSSEIDEL",
                    num_sweeps=1,
                ),
                omega_dissipation_rate_solver=SmoothSolver(
                    type="SMOOTH",
                    absolute_tolerance=1.0E-8,
                    relative_tolerance=0.01,
                    smoother="GAUSSSEIDEL",
                    num_sweeps=1,
                ),
            ),
            schemes=Schemes(
                second_order_convection=False,
            ),
            stabilization=Stabilization(
                field_limits=FieldLimits(
                    lower_temperature_bound=DimensionalTemperature(
                        value=1,
                        unit="K",
                    ),
                    upper_temperature_bound=DimensionalTemperature(
                        value=50000,
                        unit="K",
                    ),
                ),
            ),
        ),
        simulation_control=FluidSimulationControl(
            end_time=DimensionalTime(
                value=2000,
                unit="s",
            ),
            delta_t=DimensionalTime(
                value=1,
                unit="s",
            ),
            write_control=TimeStepWriteControl(
                type="TIME_STEP",
                write_interval=2000,
            ),
            num_processors=32,
            max_run_time=DimensionalTime(
                value=20000,
                unit="s",
            ),
            decompose_algorithm=ScotchDecomposeAlgorithm(
                type="SCOTCH",
            ),
        ),
        result_control=FluidResultControls(
            surface_data=[
                AreaAverageResultControl(
                    type="AREA_AVERAGE",
                    name="Top Outlet",
                    write_control=TimeStepWriteControl(
                        type="TIME_STEP",
                        write_interval=1,
                    ),
                    topological_reference=TopologicalReference(
                        entities=[
                            "B3_TE4396",
                        ],
                        sets=[],
                    ),
                ),
                AreaAverageResultControl(
                    type="AREA_AVERAGE",
                    name="Bottom Outlet",
                    write_control=TimeStepWriteControl(
                        type="TIME_STEP",
                        write_interval=1,
                    ),
                    topological_reference=TopologicalReference(
                        entities=[
                            "B11_TE4396",
                        ],
                        sets=[],
                    ),
                ),
                AreaAverageResultControl(
                    type="AREA_AVERAGE",
                    name="Top Inlet",
                    write_control=TimeStepWriteControl(
                        type="TIME_STEP",
                        write_interval=1,
                    ),
                    topological_reference=TopologicalReference(
                        entities=[
                            "B3_TE1883",
                        ],
                        sets=[],
                    ),
                ),
                AreaAverageResultControl(
                    type="AREA_AVERAGE",
                    name="Bottom Inlet",
                    write_control=TimeStepWriteControl(
                        type="TIME_STEP",
                        write_interval=1,
                    ),
                    topological_reference=TopologicalReference(
                        entities=[
                            "B11_TE1883",
                        ],
                        sets=[],
                    ),
                ),
                AreaAverageResultControl(
                    type="AREA_AVERAGE",
                    name="Battery Block 1",
                    write_control=TimeStepWriteControl(
                        type="TIME_STEP",
                        write_interval=1,
                    ),
                    topological_reference=TopologicalReference(
                        entities=[
                            "B4_TE17378",
                            "B16_TE16303",
                        ],
                        sets=[],
                    ),
                ),
                AreaAverageResultControl(
                    type="AREA_AVERAGE",
                    name="Battery Block 2",
                    write_control=TimeStepWriteControl(
                        type="TIME_STEP",
                        write_interval=1,
                    ),
                    topological_reference=TopologicalReference(
                        entities=[
                            "B4_TE17384",
                            "B16_TE16329",
                        ],
                        sets=[],
                    ),
                ),
                AreaAverageResultControl(
                    type="AREA_AVERAGE",
                    name="Battery Block 3",
                    write_control=TimeStepWriteControl(
                        type="TIME_STEP",
                        write_interval=1,
                    ),
                    topological_reference=TopologicalReference(
                        entities=[
                            "B4_TE981",
                            "B16_TE16249",
                        ],
                        sets=[],
                    ),
                ),
                AreaAverageResultControl(
                    type="AREA_AVERAGE",
                    name="Battery Block 4",
                    write_control=TimeStepWriteControl(
                        type="TIME_STEP",
                        write_interval=1,
                    ),
                    topological_reference=TopologicalReference(
                        entities=[
                            "B4_TE17372",
                            "B16_TE16223",
                        ],
                        sets=[],
                    ),
                ),
                AreaAverageResultControl(
                    type="AREA_AVERAGE",
                    name="Battery Block 5",
                    write_control=TimeStepWriteControl(
                        type="TIME_STEP",
                        write_interval=1,
                    ),
                    topological_reference=TopologicalReference(
                        entities=[
                            "B4_TE17366",
                            "B16_TE981",
                        ],
                        sets=[],
                    ),
                ),
                AreaAverageResultControl(
                    type="AREA_AVERAGE",
                    name="Battery Block 6",
                    write_control=TimeStepWriteControl(
                        type="TIME_STEP",
                        write_interval=1,
                    ),
                    topological_reference=TopologicalReference(
                        entities=[
                            "B4_TE17359",
                            "B16_TE16276",
                        ],
                        sets=[],
                    ),
                ),
            ],
            probe_points=[],
            field_calculations=[],
        ),
        contact_handling_mode="AUTO",
    )



    # Create simulation spec first. Then either use this to create a new physics-based
    # mesh, or use an existing one
    simulation_name = f"{power} per cell - Flow Rate {flowrate}"
    simulation_spec = SimulationSpec(name=simulation_name, geometry_id=geometry_id, model=model)
    simulation_id = simulation_api.create_simulation(project_id, simulation_spec).simulation_id

    # If you want to create a new mesh then add code above this line
    # For now: Get all existing meshes of the simulation project and find the existing mesh to be used
    mesh_operationS = mesh_operation_api.get_mesh_operations(project_id)


    # Take all existing meshes and find the one with the provided name at the start of the script and gets its ID.
    all_meshes = meshes_api.get_meshes(project_id)
    for m in all_meshes.embedded:
        if m.name == mesh_name:
            mesh_id = m.mesh_id
        
    # Get the simulation spec and update it with mesh_id from the previous mesh operation
    simulation_spec = simulation_api.get_simulation(project_id, simulation_id)
    simulation_spec.mesh_id = mesh_id
    simulation_api.update_simulation(project_id, simulation_id, simulation_spec)


    # Create simulation run
    simulation_run = SimulationRun(name="Run 1")
    simulation_run = simulation_run_api.create_simulation_run(project_id, simulation_id, simulation_run)
    run_id = simulation_run.run_id
    print(f"runId: {run_id}")

    # Read simulation run and update with the deserialized model
    simulation_run = simulation_run_api.get_simulation_run(project_id, simulation_id, run_id)
    simulation_run_api.update_simulation_run(project_id, simulation_id, run_id, simulation_run)

    # Start simulation run and wait until it's finished
    max_runtime = 36000
    simulation_run_api.start_simulation_run(project_id, simulation_id, run_id)
    simulation_run = simulation_run_api.get_simulation_run(project_id, simulation_id, run_id)
    simulation_run_start = time.time()
    while simulation_run.status not in ("FINISHED", "CANCELED", "FAILED"):
        if time.time() > simulation_run_start + max_runtime:
            raise TimeoutError()
        time.sleep(30)
        simulation_run = simulation_run_api.get_simulation_run(project_id, simulation_id, run_id)
        print(f"Simulation run status: {simulation_run.status} - {simulation_run.progress}")


    return {
        "result": "Hello World!",
    }