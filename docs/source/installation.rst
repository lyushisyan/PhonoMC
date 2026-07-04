Installation
============

Clone the repository
--------------------

.. code-block:: bash

   git clone https://github.com/lyushisyan/PhonoMC.git
   cd PhonoMC

Configure and build
-------------------

For a release build with OpenMP:

.. code-block:: bash

   cmake -S . -B build \
     -DCMAKE_BUILD_TYPE=Release \
     -DPHONOMC_ENABLE_OPENMP=ON \
     -DBUILD_TESTING=ON
   cmake --build build -j

The executable is written to ``build/PhonoMC``.

Run the tests
-------------

.. code-block:: bash

   ctest --test-dir build --output-on-failure

The test suite checks strict configuration parsing, boundary assignment and
periodic pairing, the provided Cross/FinFET geometry inputs, and strict POSCAR
loading.

Build without OpenMP
--------------------

.. code-block:: bash

   cmake -S . -B build -DPHONOMC_ENABLE_OPENMP=OFF
   cmake --build build -j

Build the documentation
-----------------------

.. code-block:: bash

   python3 -m pip install -r docs/requirements.txt
   make -C docs html

Open ``docs/build/html/index.html`` in a browser.
