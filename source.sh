SCRIPT="${BASH_SOURCE[${#BASH_SOURCE[@]} - 1]}"
export RNASCP_DIR="$(dirname "${SCRIPT}")"

chmod +x $RNASCP_DIR/bin/rnascp
export PYTHONPATH=$PYTHONPATH:$RNASCP_DIR/lib/python/
export PATH=$PATH:$RNASCP_DIR/bin/
export MROPATH=$RNASCP_DIR/mro/
