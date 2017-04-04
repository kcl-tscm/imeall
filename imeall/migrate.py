import os
from playhouse.migrate import *

#add_column sigma_csl to database:
DATABASE   = os.environ['GBDATABASE']
database  = SqliteDatabase(DATABASE)
database.connect()
migrator  = SqliteMigrator(database)
sigma_csl = IntegerField(default=0)
migrate(migrator.add_column('grainboundary', 'sigma_csl', sigma_csl))

