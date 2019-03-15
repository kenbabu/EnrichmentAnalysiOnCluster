import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re
import sys
import logging

logging.basicConfig( format = '%(asctime)s  %(levelname)-10s %(processName)s  %(name)s %(message)s',
                                filename='hpo_enrichment.log',
                                level = logging.INFO)

pd.set_option('display.max_colwidth', -1)



gp.__version__

names = gp.get_library_name()


def calculate_enrichment(genes,disease, bg_set='Human_Phenotype_Ontology'):
    try:
        logging.info("Start processing {}".format(disease))
        enr = gp.enrichr(gene_list= genes,
                    description='HPO_Enrichment',
                    gene_sets = bg_set,
                    outdir='HPO_HIGH_CLEANED/{}'.format(disease),
                    no_plot= True,
                    cutoff=0.5)
        logging.info("Finish processing {}".format(disease))
    except ConnectionError:
        return ("Connection Failed")
#    except Exception:
#        logging.error("Disease [{}] not processed!".format(disease))
#        pass 



def main():
    disgenet_dir = sys.argv[1]
    # disgenet_dir = '/home/user/Documents/ExternalHardDriveData/2018LITTLEMAN_EXT/2018Data'
    logging.info("Started...")
    disgenetfile = os.path.join(disgenet_dir,'Disgenet/all_gene_disease_associations.tsv')
    df_dg = pd.read_csv(disgenetfile, sep='\t', encoding = "ISO-8859-1", error_bad_lines=False)
    df_dg.head()
    groups = df_dg.groupby('diseaseId')['geneSymbol']

    dict_disgene = {k:list(v) for k,v in groups}
    df = pd.DataFrame()
    df['disease'] = dict_disgene.keys()
    df['gene'] = dict_disgene.values()
    df['num_genes'] = df['gene'].apply(lambda x: len(x))

    df_high = df.loc[(df['num_genes']>9) & (df['num_genes']<=200)]
    print(df_high.shape)
   
    ERROR_IDX = df_high.index[df_high['disease'] == 'C0280220'].tolist()[0]
    
    df_ = df_high.iloc[ERROR_IDX:]
    #
#    df_high['enrichment'] = df_high[['gene', 'disease']].apply(lambda x: calculate_enrichment(*x), axis=1 )
    df_['enrichment'] = df_[['gene', 'disease']].apply(lambda x: calculate_enrichment(*x), axis=1 )

if __name__ == '__main__':
    main()
# df = pd.DataFrame.from_dict(dict_disgene, orient='index')
