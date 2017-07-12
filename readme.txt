Getting Started

    This directory contains 4 files used in generating millions of lightcurves for MicroLensingGenerator Events. 
    The project takes input parameters from the DES Y3 data and returns a list of microlensing events for a given star object. 
    These files are:
        MicroLensingGenerator
        fake_plots
        getData
        parameters 

Prerequisites

    User must import the following into the python environment:
        import MicroLensingGenerator; import fake_plots; import getData; import matplotlib.pyplot as plt; import numpy as np; plt.ion(); import parameters
    User then imports the data from the DES survey:
        data = getData.getData(); mjd = data.get_MJD(); para = parameters.loop();
    User can then generate the plots and events:
        fake_plots.nike(mjd); obj = para.star_object(mjd)
    
    

    What things you need to install the software and how to install them

Give examples
    Installing
        A step by step series of examples that tell you have to get a development env running
    
    Say what the step will be
        Give the example
    And repeat until finished

    End with an example of getting some data out of the system or using it for a little demo

Running the tests
    Explain how to run the automated tests for this system
    Break down into end to end tests
    Explain what these tests test and why
        Give an example
            And coding style tests

    Explain what these tests test and why
        Give an example

Deployment

    Add additional notes about how to deploy this on a live system

Built With

    Dropwizard - The web framework used
    Maven - Dependency Management
    ROME - Used to generate RSS Feeds

Contributing

    Please read CONTRIBUTING.md for details on our code of conduct, and the process for submitting pull requests to us.

Versioning

    We use SemVer for versioning. For the versions available, see the tags on this repository.

Authors

    Billie Thompson - Initial work - PurpleBooth
    See also the list of contributors who participated in this project.

License

    This project is licensed under the MIT License - see the LICENSE.md file for details

Acknowledgments

    Hat tip to anyone who's code was used
    Inspiration
    etc
