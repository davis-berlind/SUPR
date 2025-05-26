#### e5_sim.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID.$TASK_ID
#$ -j y
## Job resource details:
#$ -l h_rt=3:00:00,h_data=4G
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m a
#$ -t 1-5000:1

# echo job info on joblog:
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
module load R/4.2.2

## Run code 
echo '/usr/bin/time -v Rscript ~/MICH/1e5_simulation_batch_job.R > ~/MICH/output/e5_output.$SGE_TASK_ID 2>&1' 
/usr/bin/time -v Rscript ~/MICH/1e5_simulation_batch_job.R > ~/MICH/output/e5_output.$SGE_TASK_ID 2>&1

# echo job info on joblog:
echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `date `
echo " "
#### e5_sim.sh STOP ####
