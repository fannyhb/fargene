class MetaData(object):

    def __init__(self, full_length, sens_score_file,spec_score_file,max_fpr,tmpdir,resultsdir,hmmdir):
        self.full_length = full_length
        self.sens_score_file = [sens_score_file]
        self.spec_score_file = [spec_score_file]
        self.max_fpr = max_fpr
        self.tmpdir = tmpdir
        self.resultsdir = resultsdir
        self.hmmdir = hmmdir


    def has_specificity(self):
        return self.spec_score_file != None

    def has_sensitivity(self):
        return self.sens_score_file != None
    
    def get_max_fpr(self):
        if self.max_fpr == None:
            if self.full_length:
                self.max_fpr = 0.00
            else:
                self.max_fpr = 0.05
        return self.max_fpr



