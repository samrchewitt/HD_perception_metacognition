# HD_perception_metacognition
open-access article: https://neurologyopen.bmj.com/content/4/1/e000268

Instructions for running and installing HDDM v0:8 via docker:
https://hub.docker.com/r/hcp4715/hddm

You'll need to clone this repository and launch it in docker by running command like this in terminal (Windows) 

docker run -it --rm --cpus=4 -v /c/Users/YOURUSERNAME/PATH_TO_REPO/HD_perception_metacognition/analysis/hddm_perception:/home/jovyan/hddm -p 8888:8888 hcp4715/hddm:0.8 jupyter notebook

It's also possible to run HDDM in google colab, this is the easiest way to load our models and data. However, if you want to estimate your own models, you will need Docker as Colab will timeout before models have been estimated. 

HDDM package (python) required for perceptual decision-making model : http://ski.clps.brown.edu/hddm_docs/

Metacognition model:
Hmetad package required for metacognition model implementation: https://github.com/metacoglab/HMeta-d/tree/master/Matlab

Task: 
Meta-dots task used in this study: https://github.com/metacoglab/meta_dots

Notes on dependencies:
I use plotting functions from https://github.com/anne-urai/Tools in behavioural and metacognition model scripts

Fitting models in Hmeta-D depends on JAGS 3.4.0 - see https://github.com/metacoglab/HMeta-d/tree/master/Matlab for instructions

Happy to be contacted if you have questions: s.hewitt.17@ucl.ac.uk
