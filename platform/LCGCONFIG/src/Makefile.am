includedir=$(prefix)/include/Photos

lib_LTLIBRARIES = libPhotospp.la libPhotosppHEPEVT.la

include_HEADERS =                            \
photos-C/photosC.h                           \
photos-C/pairs.h                             \
photos-C/forZ-MEc.h                          \
photos-C/forW-MEc.h                          \
photos-C/HEPEVT_struct.h                     \
photosCInterfaces/PhotosParticle.h           \
photosCInterfaces/Photos.h                   \
photosCInterfaces/PhotosEvent.h              \
photosCInterfaces/PhotosBranch.h             \
eventRecordInterfaces/PhotosHEPEVTParticle.h \
eventRecordInterfaces/PhotosHEPEVTEvent.h    \
utilities/PhotosUtilities.h                  \
utilities/PhotosRandom.h                     \
utilities/PhotosDebugRandom.h                \
utilities/Log.h

libPhotospp_la_SOURCES =             \
photos-C/photosC.cxx                 \
photos-C/forW-MEc.cxx                \
photos-C/forZ-MEc.cxx                \
photos-C/pairs.cxx                   \
photos-C/HEPEVT_struct.cxx           \
photosCInterfaces/Photos.cxx         \
photosCInterfaces/PhotosParticle.cxx \
photosCInterfaces/PhotosBranch.cxx   \
photosCInterfaces/PhotosEvent.cxx    \
utilities/Log.cxx                    \
utilities/PhotosUtilities.cxx        \
utilities/PhotosRandom.cxx
################################################################################

libPhotosppHEPEVT_la_SOURCES =                  \
eventRecordInterfaces/PhotosHEPEVTEvent.cxx     \
eventRecordInterfaces/PhotosHEPEVTParticle.cxx

if HAS_HEPMC
  lib_LTLIBRARIES += libPhotosppHepMC.la 

  libPhotosppHepMC_la_SOURCES = \
  eventRecordInterfaces/PhotosHepMCEvent.cxx    \
  eventRecordInterfaces/PhotosHepMCParticle.cxx

  include_HEADERS += \
  eventRecordInterfaces/PhotosHepMCParticle.h  \
  eventRecordInterfaces/PhotosHepMCEvent.h
  
  libPhotosppHepMC_la_LIBADD = libPhotospp.la
endif

INCLUDES = -I$(top_srcdir)/src/photos-C \
           -I$(top_srcdir)/src/photosCInterfaces \
           -I$(top_srcdir)/src/utilities \
           -I$(HEPMC_DIR)/include

libPhotosppHEPEVT_la_LIBADD = libPhotospp.la
