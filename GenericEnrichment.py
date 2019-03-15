import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re
import sys
import logging
import argparse

logging.basicConfig( format = '%(asctime)s  %(levelname)-10s %(processName)s  %(name)s %(message)s',
                                filename='generic_enrichment.log',
                                level = logging.INFO)

pd.set_option('display.max_colwidth', -1)



# gp.__version__

# names = gp.get_library_name()

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-dd", "--data_dir", help="Data directory for disease data",
                            default='/Users/Ken/2018Data',
                            # default="/home/user/Documents/ExternalHardDriveData/2018LITTLEMAN_EXT/2018Data",
                            required=False)
    parser.add_argument("-bs", "--background_set", type=int,
                            choices=[0, 1, 2, 3],

                            help=""" Mapping dictionary {0:Go_Biological_Process,
                                  1:Human_Phenotype_Ontology, 2:KEGG_2016,
                                  3: WikiPathways_2016
                                  } """
                                 )

    args = parser.parse_args()
    return args



def calculate_enrichment(genes,disease, bg_set=None):
    try:
        logging.info("Start processing {}".format(disease))
        enr = gp.enrichr(gene_list= genes,
                    description='{}_Enrichment'.format(bg_set),
                    gene_sets = bg_set,
                    outdir='{}/{}'.format(bg_set, disease),
                    no_plot= True,
                    cutoff=0.5)
        logging.info("Finish processing {}".format(disease))
    except ConnectionError:
        return ("Connection Failed")
#    except Exception:
#        logging.error("Disease [{}] not processed!".format(disease))
#        pass



def main():
    MENU = {0:"Go_Biological_Process",1:"Human_Phenotype_Ontology",
            2:"KEGG_2016",3: "WikiPathways_2016"
            }
    args = parse_arguments()
    background_set = MENU[args.background_set]
    disgenet_dir = args.data_dir

    # disgenet_dir = '/home/user/Documents/ExternalHardDriveData/2018LITTLEMAN_EXT/2018Data'
    logging.info("Started...")
    logging.info("Analysis background set: {}".format(background_set))

    disgenetfile = os.path.join(disgenet_dir,'Disgenet/all_gene_disease_associations.tsv')

    df_dg = pd.read_csv(disgenetfile, sep='\t', encoding = "ISO-8859-1", error_bad_lines=False)
    df_dg.head(2)
    groups = df_dg.groupby('diseaseId')['geneSymbol']

    dict_disgene = {k:list(v) for k,v in groups}

    df = pd.DataFrame()
    df['disease'] = dict_disgene.keys()
    df['gene'] = dict_disgene.values()
    df['num_genes'] = df['gene'].apply(lambda x: len(x))
    df_high = df.loc[(df['num_genes']>9) & (df['num_genes']<=200)]
    print(df_high.shape)
    df_ = df.iloc[:10]
    # ERROR_IDX = df_high.index[df_high['disease'] == 'C0280220'].tolist()[0]
    #
    # df_ = df_high.iloc[ERROR_IDX:]
    #
#    df_high['enrichment'] = df_high[['gene', 'disease']].apply(lambda x: calculate_enrichment(*x), axis=1 )
    df_['enrichment'] = df_[['gene', 'disease']].apply(lambda x: calculate_enrichment(*x, background_set), axis=1 )

if __name__ == '__main__':
    main()
# df = pd.DataFrame.from_dict(dict_disgene, orient='index')
