## Simple instruction for running this software

This is the folder for running a demo program for a dataset of human bodies. If you are
a MATLAB user in Windows 10, you can run the demo program by simply typing following commands:

- trainingBody (+Enter)
- runDemoBody  (+Enter)

after changing your working directory to this directory in the MATLAB command window.
The training will be executed in a leave-one-out manner.

This demo program registers the average shape and a human body included in the dataset.
If you would like to try another body shape, modify the variable `BODY_SHAPE_ID` in
`runDemoBody.m` into one of the numbers from 2 to 71.

For more information such as optional parameters, see 'README.md' above this directory,
please.
