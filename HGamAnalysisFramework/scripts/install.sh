#!/bin/bash

# Get the tag which should be checked out
tag=${1}
shift

# Get the AnalysisBase version which should be used
base=${1}
shift

# Make the proper directory structure
echo "Making packages/ directory"
mkdir packages

echo "Making run/ directory"
mkdir run # Figure out the HGamAnalysisFramework version to use
hgam_tag="HGamAnalysisFramework-$tag"
hgam_dir="${hgam_tag%%-*}/tags/${hgam_tag}"
[ "$tag" == "trunk" ]                       \
  && hgam_tag="HGamAnalysisFramework"       \
  && hgam_dir="HGamAnalysisFramework/trunk"

# Get the environment setup script
svn export svn+ssh://svn.cern.ch/reps/atlasoff/PhysicsAnalysis/HiggsPhys/Run2/HGamma/xAOD/$hgam_dir/scripts/setup_env.sh
sed -i "s/HGAM_BASE_VERSION/$base/g" setup_env.sh

# Put in the correct contact information (used by rc tag_package)
usertag="`phonebook $USER | awk '{print $2, $1}'` <`phonebook $USER --all | grep E-mail | grep -v External | awk '{print $2}'`>"
sed -i "s/HGAM_USER_TAG/$usertag/g" setup_env.sh

# Setup the environment
source setup_env.sh

cd packages/

# Get the proper HGamAnalysisFramework package
rc checkout_pkg atlasoff/PhysicsAnalysis/HiggsPhys/Run2/HGamma/xAOD/${hgam_dir}

# Get any other specified packages
for pkg in $@; do
  pkg_tag="$pkg"
  pkg_dir="${pkg_tag%%-*}/tags/${pkg_tag}"
  [ "$pkg" == "${pkg%%-*}" ] \
    && pkg_tag="$pkg"        \
    && pkg_dir="$pkg/trunk"

  rc checkout_pkg atlasoff/PhysicsAnalysis/HiggsPhys/Run2/HGamma/xAOD/${pkg_dir}
done

# Setup, find, and compile everything
./HGamAnalysisFramework/scripts/setupRelease

rc find_packages
rc compile

cd ../
