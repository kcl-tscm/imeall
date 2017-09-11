from imeall.models import GBQuery

query = GBQuery()
#modify gb_type to tilt or twist to checkout the appropriate boundary type.
query.copy_gb_dirtree(material='alphaFe', or_axis='0,0,1', 
                      pots=['PotBH.xml'], target_dir='./checkout', gb_type='tilt')

