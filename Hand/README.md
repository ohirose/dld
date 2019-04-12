## Simple instruction for running this software

This is the directory for running a demo program for the a hand dataset.
If you are a MATLAB user in Windows 10, you can run the demo program by simply typing

- trainingHand (+Enter)
- runDemoHand  (+Enter)

after changing your working directory to this directory in the MATLAB command window.
The training will be executed in a leave-one-out manner.

This demo program registers the average shape and a human handincluded in the hand
dataset. If you would like to try another hand shape, modify the variable `HAND_ID`
in `runDemoHand.m` into another number in {1 ... 40}.

For more information such as optional parameters, see `README.md` above this directory,
please.
