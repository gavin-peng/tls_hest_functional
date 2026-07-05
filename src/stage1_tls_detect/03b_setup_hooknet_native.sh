#!/bin/bash
# Set up a native Python 3.11 environment for HookNet-TLS inference (no Docker).
#
# Prerequisites (run manually first : require sudo):
#   sudo add-apt-repository ppa:deadsnakes/ppa
#   sudo apt-get update && sudo apt-get install -y python3.11 python3.11-venv python3.11-dev
#
#   # Install ASAP 2.1
#   curl --remote-name --location \
#     "https://github.com/computationalpathologygroup/ASAP/releases/download/ASAP-2.1-(Nightly)/ASAP-2.1-Ubuntu2204.deb"
#   sudo dpkg --install ASAP-2.1-Ubuntu2204.deb || sudo apt-get -f install --fix-missing --assume-yes
#   sudo ldconfig -v
#   rm ASAP-2.1-Ubuntu2204.deb
#
# Usage:
#   bash 03b_setup_hooknet_native.sh

set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../../models/pathology-hooknet-tls" && pwd)"
VENV="$SCRIPT_DIR/../../.venv_hooknet"

echo "=== Creating Python 3.11 venv at $VENV ==="
python3.11 -m venv "$VENV"

echo "=== Installing TensorFlow 2.15 + CUDA (bundled) ==="
"$VENV/bin/pip" install --upgrade pip
"$VENV/bin/pip" install "tensorflow[and-cuda]==2.15.*"

echo "=== Installing wholeslidedata + hooknet ==="
"$VENV/bin/pip" install git+https://github.com/DIAGNijmegen/pathology-whole-slide-data.git@78dd18b
"$VENV/bin/pip" install git+https://github.com/DIAGNijmegen/pathology-hooknet.git@master

echo "=== Installing remaining deps ==="
"$VENV/bin/pip" install "numpy>=1.23.5,<2.0" "scipy>=1.9.0" "shapely>=1.8.4,<2.0" tqdm pyvips

echo "=== Registering ASAP and hooknettls on Python path ==="
SITE=$("$VENV/bin/python3" -c "import site; print(site.getsitepackages()[0])")
echo '/opt/ASAP/bin' > "$SITE/asap.pth"
echo "$REPO_DIR" > "$SITE/hooknettls.pth"

echo "=== Patching hooknet private TF API ==="
find "$VENV/lib" -name "tensorflowmodel.py" \
    -exec sed -i 's|from tensorflow.python.framework.ops import Tensor|import tensorflow as tf; Tensor = tf.Tensor|' {} \;

echo ""
echo "=== Verifying GPU ==="
"$VENV/bin/python3" -c "
import tensorflow as tf
gpus = tf.config.list_physical_devices('GPU')
print(f'GPUs found: {gpus}')
"

echo ""
echo "Done. Run inference with: bash 03c_run_hooknet_native.sh NCBI688"
