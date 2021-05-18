==================
potentials-library
==================

Introduction
------------

This repository provides a copy of the records on potentials.nist.gov for use
with the potentials Python package.  This is meant as a convenience for users
who wish to work with all potentials/crystals/etc hosted in the database.

Notes:

- Downloading from this repository is faster than downloading records through
  the potentials package.  As a git repository, subsequent updates (i.e.
  pulls) are also quicker as only changes and additions need to be downloaded.

- The records are mostly in JSON format with indentations allowing them to be
  easily readable using a text editor.  The potentials package allows for the
  records to be converted to XML and/or represented in a compact way by
  copying them to a different directory if wanted/needed.

- The contained records represent a snapshop of potentials.nist.gov.  This
  should contain most or all public records in the formats listed below, but
  is not guaranteed to always be 100% up to date.

- Users are free to add their own records, clone the repository, etc.  Pull
  requests can be used to submit new records for adding to the NIST database.

Usage
-----

The directory for the potentials-library repository should correspond to the
library directory set for the potentials package:

.. code-block:: python

    import potentials
    settings = potentials.Settings()
    
    # See the current library_directory
    print(settings.library_directory)
    
    # Change the library directory (optional)
    settings.set_library_directory(PATH)

Once the records are in the library directory, the potentials.Database class
will be able to find and load them.

Contents
--------

**__schema__** contains the XSD schemas for the XML representations of
the different record types.  These schemas can be used to check and verify
if the fields inside of a record are compliant with a particular record type.

**Citation** contains the citation information associated with the
interatomic potentials stored in the database.  The citations in this
directory are saved in bibtex format, but can be converted to JSON or XML.
If the citation has a DOI, then the entry is saved using the DOI with
any forward slashes replaced with underscores to be a valid file name.  If
the citation does not have a DOI, then an ID is generated based on the
content creator(s) and year of creation.

**crystal_prototype** contains records that detail crystal prototypes, i.e.
crystal structures that have been generalized by excluding composition and
absolute dimensions.  These are used by the atomman Python package for
generating crystal structures.  The records are uniquely named through a
combination of the available associated Strukturbericht symbol, compositional
prototype and/or a common descriptor.

**dislocation** contains records that provide input parameter sets for iprPy
calculations to generate straight dislocations for common dislocation types.
The records are uniquely named based on crystal prototype, Burgers vector,
character, and slip plane family.

**free_surface** contains records that provide input parameter sets for iprPy
calculations to generate atomic systems with an ideally sliced free surface
along a crystal plane.  The records are uniquely named based on crystal
prototype and surface plane family.

**point_defect** contains records that provide input parameter sets for iprPy
calculations to generate atomic systems containing a point defect or point
defect cluster.  The records are uniquely named based on crystal prototype
and defect type.

**Potential** contains records that are used to generate the listings for
interatomic potentials found on the NIST Interatomic Potentials Repository.
These records collect citation information, usage notes, the different
available implementations, and specify the URLs where any associated
parameter files can be downloaded.  The records are uniquely named based on
citation/creation year, authors, and included interactions.

**potential_LAMMPS** contains records that allow the potentials package to
generate LAMMPS command lines for a specific LAMMPS-compatible interatomic
potential.  These records also contain the download URLs for any associated
parameter files.  The records are uniquely named based on citation/creation
year, first author, included interactions, format, and version number.

**potential_LAMMPS_KIM** contains records collecting metadata information
about openKIM models so that the potentials package can integrate installed
KIM models in with the native LAMMPS potentials.  The records are named based
on the unique KIM shortcode id.  Each record links the KIM model to one or
more Potential records and lists all known full KIM ids associated with the
different versions of the model.

**relaxed_crystal** contains records for each potential-dependent relaxed
crystal structure found by the iprPy framework for relaxations done at 0K.
The records detail the unit cell for the relaxed crystal structure which are
then useful for constructing atomic configurations using the correct lattice
constants for a given crystal structure and interatomic potential.

**stacking_fault** contains records that provide input parameter sets for
iprPy calculations to generate atomic systems containing generalized stacking
faults along a slip plane.  The records are uniquely named based on crystal
prototype and slip plane family, with "sf" added at the end to differentiate
them from the free surface records.

**related-interactions.json** lists the element and binary model interactions
for the NIST-hosted potentials and groups the related interaction models
together.  Interaction models are considered to be "related" if predictions
from the two potentials are (nearly) identical for atomic configurations
typically explored by atomistic calculations.  This characterization of
"related" is not a robust definition as there is some ambiguity if purposeful
changes to an old potential constitute an entirely new model or not.  The
remote version of this file can be obtained at
https://www.ctcms.nist.gov/potentials/site/related-interactions.json.
