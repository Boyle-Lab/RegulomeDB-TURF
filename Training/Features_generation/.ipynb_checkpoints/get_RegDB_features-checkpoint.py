import json
def get_generic_features(var_json):
    '''Function to output RegDB generic features from a json object
    Args:
        var_json(json): a json object contains RegDB features
    Returns:
        out(list of str): RegDB features ('CHIP','DNASE','PWM','FOOTPRINT','EQTL_2','PWM_matched','FOOTPRINT_matched') & ranking sores
    '''
    out = [var_json['chrom'],var_json['end']]
    feature_names = ['ChIP','DNase','PWM','Footprint','QTL','PWM_matched','Footprint_matched']
    out_features = [int(var_json['features'][feature]) for feature in feature_names] #convert to binary values
    out.extend(out_features)
    out.append(var_json['score']['ranking'])
    return out

###Organ names from new mapping file & consider multiple organs split by ',' in mapping file###
###Function to get RegulomeDB organSp features DNase & FOOTPRINT from TRIMMED COMPACT json object###
def get_organ_sp_features(var_json,organs_query,biosample_organ_mapping = 'biosample_organ.txt'):
    '''Function to output RegDB organ specific features from a json object of a variant
    Args:
        var_json(json): a json object contains RegDB features & peaks information
        organs_query (list of str): biosample organ name(s) (e.g. ['bodily fluid','blood','brain'])
        biosample_organ_mapping (str): file containing mapping between biosample name to organ name
    Returns:
        feature_values (dict of list): key: organ names in input query organ list; values: 'chrom','end','DNASE','FOOTPRINT'
    '''
    from collections import defaultdict
    import re
    biosample_organ_dict = defaultdict(str)
    #get mapping from biosample to organ
    with open(biosample_organ_mapping,'r') as f:
        for line in f:
            if len(line.strip().split('\t')) < 2: #biosample does not have corresponding organ name
                next
            else:
                biosample,organ_name = line.strip().split('\t')
                biosample_organ_dict[biosample] = organ_name #organ_name could contain multiple organs, split by ','
    ####manually add 'H1-hESC' to biosample dict (named as 'H1' in mapping file, but 'H1-hESC' in json object from .out)###
    biosample_organ_dict['H1-hESC'] = 'embryo'
    
    feature_values = defaultdict(list) #output dict
    method_names = ['DNase-seq','Footprints']
    for organ_query in organs_query:
        results = defaultdict(int) #keys: method_names; default value is 0
        for exp in var_json['peaks']:
            if exp['method'] in method_names:
                hit_organs = biosample_organ_dict[exp['biosample_term_name']].split(',') #if no mapping to organ, then return ['']; any non-empty string not in [''];
                if organ_query in hit_organs:
                    results[exp['method']] = 1
        out = [var_json['chrom'],var_json['end']]
        out.append(results['DNase-seq'])
        out.append(results['Footprints'])
        feature_values[organ_query] = out
    return feature_values
