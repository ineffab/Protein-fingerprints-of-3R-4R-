import sys

_script = 'GO_Elite.py'
_appVersion = '1.2.6'
_appDescription = "GO-Elite (http://genmapp.org/go_elite) is a software tool designed to identify a minimal non-redundant set "
_appDescription +="of Gene Ontology (GO) biological terms or pathways to describe a particular set of genes. In addition, "
_appDescription +="alternate ontologies (e.g., Disease Ontology), gene sets and metabolomics data can also be used as input."
_authorName = 'Nathan Salomonis'
_authorEmail = 'nsalomonis@gmail.com'
_authorURL = 'http://genmapp.org/go_elite'

excludes = [] #"numpy","scipy",
includes = ["suds","xml.sax.drivers2.drv_pyexpat","suds","numpy","pkg_resources","distutils"]
""" xml.sax.drivers2.drv_pyexpat is an XML parser needed by suds that py2app fails to include. Identified by looking at the line: parser_list+self.parsers in
/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/PyXML-0.8.4-py2.7-macosx-10.6-intel.egg/_xmlplus/sax/saxexts.py
check the py2app print out to see where this file is in the future"""

if sys.platform.startswith("darwin"):
        from distutils.core import setup
        import py2app
        options = {"py2app":
                    {"excludes": excludes,
                     "includes": includes,
                     "compressed":True,
                    "iconfile": "goelite.icns"}
        }
        setup(name="GO_Elite",
                        app=[_script],
                        version=_appVersion,
                        description=_appDescription,
                        author=_authorName,
                        author_email=_authorEmail,
                        url=_authorURL,
                        options=options,
                        #data_files=data_files,
                        setup_requires=["py2app"]
        )
