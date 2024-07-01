# Short Description
This repository hosts the API python scripts used for the joint SimScale & Valispace webinar.

## Available files and workflow

There are two example script files provided

You can use the [01_SimScale-Valispace_batteryPack_runSimulation_doc.py](https://github.com/SimScaleGmbH/simscale-valispace-api/blob/master/01_SimScale-Valispace_batteryPack_runSimulation_doc.py) to see how to:

* Import relevant SimScale and Valispace Libraries
* How to read _Valis_ from your Valispace project and pass them to SimScale.
* How to read, edit and create a new specs for a Conjugate Heat Transfer analysis in SimScale
* How to start a new simulation in SimScale via the API
  
Once the simulation is finished, you can utilize the [02_SimScale-Valispace_batteryPack_postProcessingOnly_doc.py](https://github.com/SimScaleGmbH/simscale-valispace-api/blob/master/02_SimScale-Valispace_batteryPack_postProcessingOnly_doc.py) to learn how to:

* Import the relevant post-processing SimScale libraries
* Retrieve simulation results like plots, 3D fields, screeshots etc
* Feed the simulation results back to Valispace by updating the _Valis_
* Update an Analysis Report in Valispace w/ new screenshots

A sample project where those scripts can be applied is the following
https://www.simscale.com/projects/itsavlidis_api/webinar-_battery_pack_cooling_-_design_study_-_valispace/

>[!TIP]
>The scripts above can be used as a quick templates. If you are particularly interested in creating a scirpted SimScale-Valispace workflow for your own applications please contact us.

## SimScale-Valispace integration

You can learn about the **SimScale/Valispace integration** by watching our webinar on-demand:

_A Simple AI & CFD Workflow for EV Battery Pack Cooling Desing_
https://www.simscale.com/webinars-workshops/a-simple-ai-cfd-workflow-for-ev-battery-pack-cooling-design/

or visting our integration page
https://www.simscale.com/product/integrations-partners/valispace/


## Benefits
* See how design changes influence your product behavior.
* Updating any Valispace field seamlessly pushes a new simulation to SimScale and results are rapidly fed back, taking advantage of SimScale’s cloud computing.
* Build automated simulation reports and go/no-go requirements matrices to understand your system’s behavior.
* Enable your design teams to gain a rapid understanding of design changes by performing high-fidelity SimScale simulations in the background with a push of a button.


## Requirements
* SimScale Professional license
* SimScale API key
* Valispace Pro or Enterprise license
