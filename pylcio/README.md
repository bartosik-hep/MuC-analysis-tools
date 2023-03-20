# PyLCIO drivers

Different `Driver` implementations in the `drivers` folder take an `LCIO` event as input and fill a `ROOT::TTree` as output.
Plots can be produced from the resulting trees in a preferred way, including Jupyter Notebooks collected in [`/notebooks`](/notebooks/).

Edit `run.py` to import the driver of interest and run it over the input `*.slcio` files.

PyLCIO provides high flexibility at the expense of much slower performance compared to a compiled Marlin processor in C++.
