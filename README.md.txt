## PIBIC 2024–2025

Este repositório contém os códigos do projeto PIBIC 2024–2025.

## Referências usadas

- [Algoritmo AprilTag (versão inicial)](https://people.csail.mit.edu/kaess/apriltags/)
- [Detecção ArUco com OpenCV](https://docs.opencv.org/4.x/d5/dae/tutorial_aruco_detection.html)
- [Calibração com OpenCV em Python](https://docs.opencv.org/4.x/dc/dbb/tutorial_py_calibration.html)

Além disso, para utilizar a **remoteAPI no CoppeliaSim**, foi necessário instalar a biblioteca *apriltag* no Linux — o processo padrão de instalação copia os arquivos para `usr/local/include` e `usr/local/lib`, além de gerar um script pkg-config em `usr/local/lib/pkgconfig`.

No cliente C++ via ZMQ (dentro da pasta `programming/zmqRemoteApi/clients/cpp`), utilizei uma configuração CMake similar à oferecida por padrão. As dependências `jsoncons` e `cppzmq` são baixadas automaticamente, mas podem ser instaladas globalmente para otimizar builds futuros. O processo utilizado foi:

```bash
mkdir build
cd build
cmake -DGENERATE=OFF ..
cmake --build . --config Release
