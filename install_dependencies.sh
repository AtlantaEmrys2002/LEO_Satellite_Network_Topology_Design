# Installs package dependencies on a macOS-based system

echo "INSTALL PROJECT DEPENDENCIES"

echo "Please note that users must have pip installed on their machine."

echo "Installing dependencies... "
echo ""

# Python packages
pip install astropy certifi ephem networkx numpy nx-parallel pytest scipy sgp4 skyfield sphinx || exit 1

echo "Dependencies successfully installed"
