# Installs package dependencies on a macOS-based system

echo "INSTALL PROJECT DEPENDENCIES"

echo "Please note that users must have pip installed on their machine."

echo "Installing dependencies... "
echo ""

# Python packages
pip install astropy certifi configobj ephem matplotlib networkx numpy plotly pytest scipy sgp4 skyfield sphinx || exit 1
pip install --upgrade scipy  # Ensures scipy.sparse.csr_array installed correctly
pip install --upgrade certifi

echo "Dependencies successfully installed"
