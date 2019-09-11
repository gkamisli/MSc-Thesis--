import numpy as np
# FIND IN WHICH GENES GIVEN SNPS ARE
nNs = np.array( [vBim_pNC['snp'][i] for i in range(len(vBim_pNC))] ) # snp ids from vBim_pNC column

print( "Find gene ids of HBTRC study snps" )
from Bio import Entrez
Entrez.email = "email address"
gene_snp_records = np.array( [Entrez.read(Entrez.elink(dbfrom="snp", 
                                  id=snp_ids.replace('rs', ''), 
                                  db="gene")) for snp_ids in nNs] )

# for more than 10 gene records, we need to tackle genes that are not related to the given snps
gene_snp_uids = np.array( [ gene_record[0]['LinkSetDb'][0]['Link'][0]['Id']
                           for gene_record in gene_snp_records ] )

print( "Find Chromosome locations of HBTRC examples" )
iChromosome_nX = np.full( shape = gene_snp_uids.shape,
                          fill_value = np.nan )

# Request data for all your genes  
request = Entrez.epost( "gene",
                        id = ",".join( gene_snp_uids ) )
result = Entrez.read(request)
webEnv = result[ "WebEnv" ]
queryKey = result[ "QueryKey" ]
data = Entrez.esummary( db = "gene", 
                        webenv = webEnv, 
                       query_key = queryKey )
annotations = Entrez.read( data )

# Wrangle chromosome locations (uid: unique identifier)
sGeneID_0_nG = np.array( [ gene.attributes[ "uid" ] for gene in annotations[ "DocumentSummarySet" ][ "DocumentSummary" ] ] )
iChromosome_0_nG = np.array( [ gene[ "Chromosome" ] for gene in annotations[ "DocumentSummarySet" ][ "DocumentSummary" ] ] )
assert np.all( gene_snp_uids == sGeneID_0_nG )

# The start and ending locations of each gene on chromosome
iChrStart_0_nG = np.array( [ gene["GenomicInfo"][0]["ChrStart"] for gene in annotations[ "DocumentSummarySet" ][ "DocumentSummary" ] ] )
iChrStop_0_nG = np.array( [ gene["GenomicInfo"][0]["ChrStop"] for gene in annotations[ "DocumentSummarySet" ][ "DocumentSummary" ] ] )

print( "Location search done.")