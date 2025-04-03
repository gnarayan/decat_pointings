This GH repo contains the observing plans, logs and scripts for the DECam Alliance For Transients. 
You will find a copy of this repo on the CTIO DECam machine, observer4

Every observing night:
- Talk to the Observing Assistant if you need information on how to connect to the CTIO DECam machines
  - do this early, or you will not be observing as no one but CTIO staff can give you this info 
- Clone/update this repo
  - `git pull` if an existing machine or `git@github.com:gnarayan/decat_pointings.git` if on a new machine
- Update the observing plans if ncessary, and make sure your program's observing scripts are in the `jsons/` directory
        - if you are on observer4 and for some reason need to check the repo out again then
            - `cd /home/DECamObserver/ExposureScripts/User_scripts`
            - `git clone git@github-decat:gnarayan/decat_pointings.git`
            - NOTE THE `github-decat` alias instead of `github.com`!!!
- Observe at night - detailed instructions are on the DECAT slack workspace - contact Gautham Narayan <gsn@illinois.edu> if you need to be added 
(or ask your program PI) to the slack
- Run `qcInvPrint` after the last observation has been processed by the DECAM Data Handling System (DHS)
- Use the `upload_qcinv.sh` script to upload the .qcinv logfile
- Complete and submit the night report for the observatory staff
