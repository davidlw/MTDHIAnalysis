if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    from CRABClient.UserUtilities import config, getUsernameFromSiteDB
    config = config()

    config.General.workArea = 'Gen'
    config.General.transferOutputs = True
    config.General.transferLogs = False
    config.JobType.pluginName = 'Analysis'
    config.JobType.numCores = 4
    config.JobType.maxMemoryMB = 8000
    config.JobType.maxJobRuntimeMin = 2750
    config.Data.unitsPerJob = 1
#    config.Data.totalUnits = 10
    config.Data.splitting = 'FileBased'
#    config.Data.splitting = 'Automatic'
#    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/flowcorr/'
#    config.Data.publication = True
#    config.Data.useParent = True
    config.Data.inputDBS = 'phys03'
    config.Site.storageSite = 'T2_CH_CERN'

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    config.General.requestName = 'Hydjet_Gen_LambdaC_TOF_woETL_v6'
    config.JobType.psetName = '../test/test_woETL_lamC_cfg.py'
    config.Data.inputDataset = '/Hydjet_Quenched_MinBias_5020GeV_PhaseII/davidlw-GEN_SIM_102X_test_v2-9d87269d5e2b5ba986bce8e2a858aaf3/USER'
    config.Data.outputDatasetTag = 'Gen_LambdaC_TOF_woETL_v6'
    submit(config)

    config.General.requestName = 'Hydjet_Gen_LambdaC_TOF_wETL_v6'
    config.JobType.psetName = '../test/test_wETL_lamC_cfg.py'
    config.Data.inputDataset = '/Hydjet_Quenched_MinBias_5020GeV_PhaseII/davidlw-GEN_SIM_102X_test_v2-9d87269d5e2b5ba986bce8e2a858aaf3/USER'
    config.Data.outputDatasetTag = 'Gen_LambdaC_TOF_wETL_v6'
    submit(config)
