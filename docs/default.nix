with import <nixpkgs> { };

let
  venvDir = "./.venv";
  pythonPackages = pkgs.python310Packages;
in pkgs.mkShell rec {
  name = "impurePythonEnv";
  buildInputs = [
    pythonPackages.python
    pythonPackages.virtualenv
    pythonPackages.numpy
    pkgs.zlib
    pkgs.expat
    nodejs_20
  ];

  # This is very close to how venvShellHook is implemented, but
  # adapted to use 'virtualenv'
  shellHook = ''
    SOURCE_DATE_EPOCH=$(date +%s)

    if [ -d "${venvDir}" ]; then
      echo "Skipping venv creation, '${venvDir}' already exists"
    else
      echo "Creating new venv environment in path: '${venvDir}'"
      ${pythonPackages.python.interpreter} -m venv "${venvDir}"
    fi

    source "${venvDir}/bin/activate"

    pip install -r requirements.txt

    LD_LIBRARY_PATH="${pkgs.lib.makeLibraryPath buildInputs}:$LD_LIBRARY_PATH"
    LD_LIBRARY_PATH="${pkgs.stdenv.cc.cc.lib.outPath}/lib:$LD_LIBRARY_PATH"

    # Add the expat library path to LD_LIBRARY_PATH
    LD_LIBRARY_PATH=${pkgs.expat}/lib:$LD_LIBRARY_PATH

    # fixes libstdc++ issues and libgl.so issues
    LD_LIBRARY_PATH=${pkgs.stdenv.cc.cc.lib}/lib/:/run/opengl-driver/lib/:$LD_LIBRARY_PATH
    # fixes xcb issues :
    QT_PLUGIN_PATH=${pkgs.qt5.qtbase}/${pkgs.qt5.qtbase.qtPluginPrefix}
    PATH="$PWD/node_modules/.bin/:$PATH"
  '';
}
