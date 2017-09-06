from imeall.models import GBQuery

query = GBQuery()
query.copy_gb_dirtree(material='alphaFe', or_axis='0,0,1', pots=['PotBH.xml'], target_dir='./checkout')

