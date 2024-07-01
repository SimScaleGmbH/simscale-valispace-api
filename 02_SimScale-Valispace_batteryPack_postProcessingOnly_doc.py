from typing import Any, Dict
import requests
from io import BytesIO

import valispace

# Import SimScale required libraries
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
from simscale_sdk import AbsolutePowerSource, SmoothSolver

# Simscale postprocessing libraries
from simscale_sdk import ScalarField, ModelSettings, CuttingPlane, Vector3D, RenderMode, Filters, UserInputCameraSettings, ProjectionType, ScreenshotOutputSettings, ResolutionInfo, ScreenshotReportProperties
from simscale_sdk import ReportRequest, Part

import time
import urllib3
import csv

import numpy as np
import pandas as pd
import plotly.io as pio
import plotly.graph_objects as go

# Importing the SIMSCALE_API_KEY from user secrets defined in Settings
from .settings import simscale_key

# Importing Valispace username and password to be used in the script
from .settings import Username, Password

deployment = "simscale"
valispace = valispace.API(url=f'https://{deployment}.valispace.com', username=Username, password=Password)
projectId = 21

# This is the function used to massively upload created images to the Analysis report
def upload_images(deployment: str, imageIndex:dict, figureIndex: dict,  filePath: dict):
    
    for key in imageIndex:
        # Get the file path of the saved image
        file_path = filePath[key]

        # Read the image file
        with open(file_path, 'rb') as file:
            image_data = file.read()

        # Upload Analysis image
        url = f"https://{deployment}.valispace.com/rest/analyses/blocks/image/{figureIndex[key]}/"

        files = [
            ('file', (f'{key}.png', image_data, 'image/png'))
        ]
        headers = {
            'Authorization': valispace._session.headers['Authorization']
        }

        response = requests.request("PATCH", url, headers=headers, files=files)

    print(response.text)


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
    simulation_api = SimulationsApi(api_client)
    simulation_run_api = SimulationRunsApi(api_client)
    reports_api = ReportsApi(api_client)
    #materials_api = MaterialsApi(api_client)


    



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



    # Start Results PP

    # This is a decoupled script. We do the post-processing as a separate step from the simulation setup and run.
    # However those can be easily combined into a single script

    # For this demonstration we can pick the simulation to postprocess by its name.

    # This script only updates the OUTPUT values and the analysis Screenshots and Graphs
    # It needs to be run in combination with the "SimScale_batteryPack_runSimulation" script along with
    # some modifications to consistently update both inputs and outputs

    # All available completed simulations in the demo project
    PP_sim_1 = "Existing Design -12W per cell 0.13 kg/s"
    PP_sim_2 = "13.0W per cell - Flow Rate 0.13kg/s"
    PP_sim_3 = "13.0W per cell - Flow Rate 0.26kg/s"

    # Pick which one to post-process 
    PP_sim = PP_sim_2

    # Get simulations of the current project and find the one with the given name.
    simulations = simulation_api.get_simulations(project_id)

    for sim in simulations.embedded:
        if sim.name == PP_sim:
            simulation_id = sim.simulation_id

    simulation_runs = simulation_run_api.get_simulation_runs(project_id,simulation_id)
    run_id = simulation_runs.embedded[0].run_id


    # Caclulate run duration, core-hour consumption and cost
    run_info = simulation_run_api.get_simulation_run(project_id,simulation_id,run_id)
    run_duration = run_info.finished_at - run_info.created_at
    run_duration_minutes = round(run_duration.seconds/60)

    run_spec = simulation_run_api.get_simulation_run_spec(project_id,simulation_id,run_id)
    # Get machince cpus from simulation specs and calculate corehour consumption and total cost in euros
    run_machine_cpus = run_spec.model.simulation_control.num_processors
    run_cpuh = round(run_machine_cpus*run_duration_minutes/60,2)
    cpuh_cost = 0.1 #Cost may change based on the size of the corehour package purschased. This can be also set directly in valispace.
    total_cost = round(cpuh_cost*run_cpuh,2)

    print(PP_sim)
    print(f"This Run took {run_duration_minutes} minutes to compute and consumed {run_cpuh} corehours."
        f"The total cost of this run is {total_cost} euros.")

    # Get full results info
    results_full = simulation_run_api.get_simulation_run_results(project_id, simulation_id, run_id)


    # This function is used to retrieve the csv data for any plot-like result available in SimScale.
    def get_csv_data(category, name, quantity):
        
        results = simulation_run_api.get_simulation_run_results(project_id, simulation_id, run_id, category = category, name = name, quantity = quantity)
        results = results.embedded[0]
        results_response = api_client.rest_client.GET(url=results.download.url,
                                                    headers={api_key_header: api_key},
                                                    _preload_content=False)
        results_csv = results_response.data.decode("utf-8")

        results_csv = results_csv.splitlines()
        results_csv = csv.reader(results_csv,delimiter =",")
        results_csv = list(results_csv)
        
        return results_csv




    # Get pressure drop
    results_dp = get_csv_data(None,"Top Inlet", "p")
    dp = round(float(results_dp[-1][1]))

    print(f"The pressure drop is {dp} Pa.")



    # Get the maximum average and the minimum average temperature per block.
    # This block is added to demonstrate how to get averaged min/max temperature values from SimScale result controls.
    # For the overall min/max/ave values a 3rd party tool i.e. Paraview needs to be used at the moment if the result are
    # meant to be retrieved programmatically.

    # Currently these values i.e. overall min/max/ave can be manually retrieved through the post-processor GUI via "Statistics"
    # however this might change in the future, as we plan to add direct results controls to retrieve those values as well
    # which would be able to be accessed programmatically.


    # This block of code is not actually used - it is just to demo how-to retrieve the data from plots if needed.
    # The actual used values are hardcoded later
    def min_max_block_T(blockNumber, data):    
        min_ave = round(min(float(data[-1][1]),float(data[-1][2])) - 273.15, 2)
        max_ave = round(max(float(data[-1][1]),float(data[-1][2])) - 273.15, 2)
        
        return [blockNumber, min_ave, max_ave]


    numOfBlocks = 6
    block_Average_MinMax = [] 

    for i in range(1,numOfBlocks+1):
        block_data = get_csv_data(None,f"Battery Block {i}", "T")
        block_Average_MinMax.append(min_max_block_T(i, block_data))

    print("The min average and the max average temperature per block are:")
    for r in block_Average_MinMax:
        print(r)
    # End of block of code (not used)


    # Block of code to demonstrate how to update valis based on simulation results.

    #****** VALISPACE UPDATE VALUES **********

    #update duration time in minutes
    temp_vali_ID = 12908
    update = valispace.request('patch',f"valis/{temp_vali_ID}/", {"formula": f"{run_duration_minutes}min"})

    #update number of cpus 
    temp_vali_ID = 12909
    update = valispace.request('patch',f"valis/{temp_vali_ID}/", {"formula": f"{run_machine_cpus}"})

    #update the total pressure drop value 
    temp_vali_ID = 12912
    update = valispace.request('patch',f"valis/{temp_vali_ID}/", {"formula": f"{dp}Pa"})


    # Hardcoded min/max/ave cell values i.e. 
    # Currently these values i.e. overall min/max/ave can be manually retrieved through the post-processor GUI via "Statistics"
    # or via a 3rd Party tool (Paraview)
    # however this might change in the future, as we plan to add direct result controls to retrieve those values as well
    # which would be able to be accessed programmatically.

    if PP_sim == PP_sim_1:
        [T_min, T_max, T_ave] = [21.26, 37.8, 29.96]
    elif PP_sim == PP_sim_2:
        [T_min, T_max, T_ave] = [21.52, 39.44, 30.95]
    elif PP_sim == PP_sim_3:
        [T_min, T_max, T_ave] = [20.93, 38.36, 29.94]

    # update T_min  
    temp_vali_ID = 12878
    update = valispace.request('patch',f"valis/{temp_vali_ID}/", {"formula": f"{T_min}degC"})

    # update T_max  
    temp_vali_ID = 12876
    update = valispace.request('patch',f"valis/{temp_vali_ID}/", {"formula": f"{T_max}degC"})

    # update T_ave  
    temp_vali_ID = 12877
    update = valispace.request('patch',f"valis/{temp_vali_ID}/", {"formula": f"{T_ave}degC"})


    #****** VALISPACE UPDATE VALUES END **********


    # Image post-processing

    # Create dictionaries to host image data and then submitted to image function
    imageIndex = {}
    figureIndex = {}
    filePath = {}


    # Create Domain Residuals plot - 1st image
    pd.options.plotting.backend = "plotly"
    pio.renderers.default='browser'

    results_domRes_csv = get_csv_data("DOMAIN_PLOT", None, None)

    df = pd.DataFrame(results_domRes_csv)
    new_header = df.iloc[0].tolist()
    df = df[1:]
    df.columns = new_header
    df = df.astype(float)

    fig = df.plot("Time (s)", ['T', 'Ux', 'Uy', 'Uz', 'p'], log_y=False).update_layout(
        yaxis_title = "Domain Convergence",
        legend_title = None,
        plot_bgcolor='white',
        legend=dict(font=dict(size= 15))
        )
    fig.update_xaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey',
        title_font = dict(size=20)
    )
    fig.update_yaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey',
        title_font = dict(size=20)
    )
    #fig.show()

    # Populate dictonaries
    current_image_key = "Domain Residuals"
    imageIndex[current_image_key] = 1
    figureIndex[current_image_key] = 96
    file_path = './01_domain_plotly.png'
    filePath[current_image_key] = file_path

    # Save image
    try:
        print("Attempting to save the image...")
        pio.write_image(fig, file_path, format='png')
        print(f"Image successfully saved at {file_path}")
    except Exception as e:
        print(f"An error occurred: {e}")


    # Create Residuals Plot **** - 2nd image

    results_res_csv = get_csv_data("RESIDUALS_PLOT", None, None)

    df =[]

    df = pd.DataFrame(results_res_csv)
    new_header = df.iloc[0].tolist()
    df = df[1:]
    df.columns = new_header
    df = df.astype(float)



    fig = df.plot("Time (s)", ['T', 'Ux', 'Uy', 'Uz', 'k', 'omega', 'p_rgh'], log_y=True).update_layout(
        yaxis_title = "Residuals",
        legend_title = None,
        plot_bgcolor='white',
        legend=dict(font=dict(size= 15))
        )
    fig.update_xaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey',
        title_font = dict(size=20),
        rangemode = "tozero"
    )
    fig.update_yaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey',
        title_font = dict(size=20),
        range = [-7,0],
        tickformat = '.1e',
        tickvals=[1,0.1,0.01,0.001,0.0001,0.00001,0.000001,0.0000001]
    )

    current_image_key = "Residuals"
    imageIndex[current_image_key] = 2
    figureIndex[current_image_key] = 97
    file_path = './02_residuals_plotly.png'
    filePath[current_image_key] = file_path

    try:
        print("Attempting to save the image...")
        pio.write_image(fig, file_path, format='png')
        print(f"Image successfully saved at {file_path}")
    except Exception as e:
        print(f"An error occurred: {e}")



    #******************************

    # 3D results PP and screenshot creation

    results_3D = simulation_run_api.get_simulation_run_results(project_id, simulation_id, run_id, category = "SOLUTION")
    solution_info = results_3D.embedded[0]


    # # Simulation Report

    report_request = []
    create_report_response = []



    def reporting(project_id, report_id, screenShot_Name):
        print(f"Starting report with ID {report_id} for screenshot {screenShot_Name}")
        report_job = reports_api.start_report_job(project_id, report_id)

        report = reports_api.get_report(project_id, report_id)

        while report.status not in ("FINISHED", "CANCELED", "FAILED"):
            time.sleep(30)
            report = reports_api.get_report(project_id, report_id)

        print(f"Report finished with status {report.status}")

        if report.status == "FINISHED":
            # Download the report
            print("Downloading report result")
            report_response = api_client.rest_client.GET(
                url=report.download.url,
                headers={api_key_header: api_key},
                _preload_content=False,
            )

            file_name = f"{screenShot_Name}.{report.download.format}"
            with open(file_name, "wb") as file:
                file.write(report_response.data)
                print(f"Finished downloading report with name {file_name}")
        elif report.status == "FAILED":
            raise Exception("Report generation failed", report.failure_reason)
            
        return []




    #*** 03 Screenshot***
    screenShot_Name = "03_velocity"

    scalar_field = ScalarField(field_name="Velocity", component ="Magnitude", data_type="CELL")


    model_settings = ModelSettings(
        parts=[],
        scalar_field=scalar_field
        
    )

    cutting_plane = CuttingPlane(
        name="Velocity-plane",
        scalar_field=scalar_field,
        center=Vector3D(x=-0.04853767156600952, y=0.13077771663665771, z=0.06675655394792557),
        normal=Vector3D(x=0, y=0, z=-1),
        opacity=1,
        clipping=True,
        render_mode=RenderMode.SURFACES
    )

    filters = Filters(cutting_planes=[cutting_plane])

    camera_settings = UserInputCameraSettings(
        projection_type=ProjectionType.ORTHOGONAL,
        up=Vector3D(0, -1, 0),
        eye=Vector3D(-0.25264, 0.258597, 0.6509),
        center=Vector3D(-0.2526, 0.2585, -0.349005),
        front_plane_frustum_height=0.497228,
    )



    output_settings = ScreenshotOutputSettings(
        name="Output 1",
        format="PNG",
        resolution=ResolutionInfo(x=1440, y=1080),
        frame_index=1,
        show_legend=True,
        show_cube=False
    )
    report_properties = ScreenshotReportProperties(
        model_settings=model_settings,
        filters=filters,
        camera_settings=camera_settings,
        output_settings=output_settings,
    )
    report_request.append(ReportRequest(
        name="Report 1",
        description="Simulation report",
        result_ids=[solution_info.result_id],
        report_properties=report_properties
    ))


    create_report_response.append(reports_api.create_report(project_id, report_request[0]))
    report_id = create_report_response[0].report_id

    # Start report job


    current_image_key = "Velocity"
    imageIndex[current_image_key] = 3
    figureIndex[current_image_key] = 98
    file_path = f'./{screenShot_Name}.PNG'
    filePath[current_image_key] = file_path

    reporting(project_id, report_id, screenShot_Name)

    #*** 04 Screenshot***
    screenShot_Name = "04_pressure"

    scalar_field = ScalarField(field_name="Pressure", data_type="CELL")
    model_settings = ModelSettings(
        parts=[],
        scalar_field=scalar_field
        
    )
    cutting_plane = CuttingPlane(
        name="Pressure-plane",
        scalar_field=scalar_field,
        center=Vector3D(x=-0.04853767156600952, y=0.13077771663665771, z=0.06675655394792557),
        normal=Vector3D(x=0, y=0, z=-1),
        opacity=1,
        clipping=True,
        render_mode=RenderMode.SURFACES
    )
    filters = Filters(cutting_planes=[cutting_plane])

    camera_settings = UserInputCameraSettings(
        projection_type=ProjectionType.ORTHOGONAL,
        up=Vector3D(0, -1, 0),
        eye=Vector3D(-0.25264, 0.258597, 0.6509),
        center=Vector3D(-0.2526, 0.2585, -0.349005),
        front_plane_frustum_height=0.497228,
    )
    output_settings = ScreenshotOutputSettings(
        name="Output 1",
        format="PNG",
        resolution=ResolutionInfo(x=1440, y=1080),
        frame_index=1,
        show_legend=True,
        show_cube=False
    )
    report_properties = ScreenshotReportProperties(
        model_settings=model_settings,
        filters=filters,
        camera_settings=camera_settings,
        output_settings=output_settings,
    )
    report_request.append(ReportRequest(
        name="Report 1",
        description="Simulation report",
        result_ids=[solution_info.result_id],
        report_properties=report_properties
    ))


    create_report_response.append(reports_api.create_report(project_id, report_request[1]))
    report_id = create_report_response[1].report_id


    current_image_key = "Pressure"
    imageIndex[current_image_key] = 4
    figureIndex[current_image_key] = 99
    file_path = f'./{screenShot_Name}.PNG'
    filePath[current_image_key] = file_path

    reporting(project_id, report_id, screenShot_Name)




    #*** 05 Screenshot***
    screenShot_Name = "05_Internal Cell Temperature"

    scalar_field = ScalarField(field_name="Temperature", data_type="CELL")
    model_settings = ModelSettings(
        parts=[],
        scalar_field=scalar_field
        
    )
    cutting_plane = CuttingPlane(
        name="Temperature-plane",
        scalar_field=scalar_field,
        center=Vector3D(x=-0.06300857663154602, y=0.3658024072647095, z=0.021883537992835045),
        normal=Vector3D(x=0, y=0, z=-1),
        opacity=1,
        clipping=True,
        render_mode=RenderMode.SURFACES
    )
    filters = Filters(cutting_planes=[cutting_plane])

    camera_settings = UserInputCameraSettings(
        projection_type=ProjectionType.ORTHOGONAL,
        up=Vector3D(-0.22115817666053772, -0.19899803400039673,0.9547192454338074),
        eye=Vector3D(0.27182525396347046, 0.6647929549217224, 0.2565128803253174),
        center=Vector3D(-0.47261548042297363, 0.06683850288391113, -0.040570199489593506),
        front_plane_frustum_height=0.4550038415599161,
    )
    output_settings = ScreenshotOutputSettings(
        name="Output 1",
        format="PNG",
        resolution=ResolutionInfo(x=1440, y=1080),
        frame_index=1,
        show_legend=True,
        show_cube=False
    )
    report_properties = ScreenshotReportProperties(
        model_settings=model_settings,
        filters=filters,
        camera_settings=camera_settings,
        output_settings=output_settings,
    )
    report_request.append(ReportRequest(
        name="Report 1",
        description="Simulation report",
        result_ids=[solution_info.result_id],
        report_properties=report_properties
    ))


    create_report_response.append(reports_api.create_report(project_id, report_request[2]))
    report_id = create_report_response[2].report_id

    current_image_key = "Internal Temperature"
    imageIndex[current_image_key] = 5
    figureIndex[current_image_key] = 100
    file_path = f'./{screenShot_Name}.PNG'
    filePath[current_image_key] = file_path
    reporting(project_id, report_id, screenShot_Name)



    #*** 06 Screenshot***
    screenShot_Name = "06_Overall Battery Temperature"

    scalar_field = ScalarField(field_name="Temperature", data_type="CELL")

    hidden_parts = [Part(part_identifier="region180"),
                    Part(part_identifier="region181"),
                    Part(part_identifier="region170"),
                    Part(part_identifier="region172"),
                    Part(part_identifier="region171"),
                    Part(part_identifier="region171")]


    model_settings = ModelSettings(
        parts=hidden_parts,
        hide_selected_parts=True,
        show_volumes = False,
        scalar_field=scalar_field
        )



    camera_settings = UserInputCameraSettings(
        projection_type=ProjectionType.ORTHOGONAL,
        up=Vector3D(-0.22115817666053772, -0.19899803400039673,0.9547192454338074),
        eye=Vector3D(0.27182525396347046, 0.6647929549217224, 0.2565128803253174),
        center=Vector3D(-0.47261548042297363, 0.06683850288391113, -0.040570199489593506),
        front_plane_frustum_height=0.4550038415599161,
    )
    output_settings = ScreenshotOutputSettings(
        name="Output 1",
        format="PNG",
        resolution=ResolutionInfo(x=1440, y=1080),
        frame_index=1,
        show_legend=True,
        show_cube=False
    )
    report_properties = ScreenshotReportProperties(
        model_settings=model_settings,
        filters=None,
        camera_settings=camera_settings,
        output_settings=output_settings,
    )
    report_request.append(ReportRequest(
        name="Report 1",
        description="Simulation report",
        result_ids=[solution_info.result_id],
        report_properties=report_properties
    ))


    create_report_response.append(reports_api.create_report(project_id, report_request[3]))
    report_id = create_report_response[3].report_id

    current_image_key = "Overall temperature"
    imageIndex[current_image_key] = 6
    figureIndex[current_image_key] = 101
    file_path = f'./{screenShot_Name}.PNG'
    filePath[current_image_key] = file_path
    reporting(project_id, report_id, screenShot_Name)


    # Once all images are created use this function to upload them to the analysis report.
    upload_images(deployment, imageIndex, figureIndex, filePath)

    return {
        "result": "Hello World!",
    }