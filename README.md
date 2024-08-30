# 20Hz Data Analysis description
This is a re-analysis of a behavioral task performance where parameters of deep brain stimulation (DBS) in patients with Parkinson's disease were systematically varied. Participants performed a simple reaction time task, a flanker task, a modified version of the Go-NoGo task, and a signal-change task (variation of the stop-signal task). Participants performed the tasks while receiving their usual high-frequency stimulation to treat their motor symptoms, our experimental low-frequency 20Hz stimulation or no stimulation (OFF). Both the experimenter and participants were blinded to the stimulation condition.

We use Bayesian mixed effect models with a shifted-log normal likelihood for the reaction time analysis and a Bernoulli likelihood with a logit link function to analyze the error data. Mixed effect modeling is performed in R using the brms package.

The subfolder "Analysis" contains the scripts for the described analyses, and the folder "Figures" contains scripts to create the visualizations, mainly performed using ggplot in R.

To use scripts and load model files use package versions contained in the renv file. Bayes_factor calculations may depend on brms version or dependency used to create respective models. To use renv files:

```
# Install renv if not already installed
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

# Restore the project environment
renv::restore()
```
