version = {}
with open("...imeall/version.py") as fp:
    exec(fp.read(), version)
# ...
setup(version = version['__version__'])
