#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import h5py

def explore_hdf5(filename):
    with h5py.File(filename, "r") as f:
        def print_attrs(name, obj):
            print(f"{name:30s} shape={obj.shape} dtype={obj.dtype}")
        f.visititems(print_attrs)

if __name__ == "__main__":
    # 注意这里路径要加引号，并且是字符串
    explore_hdf5("kappa-m31311.hdf5")