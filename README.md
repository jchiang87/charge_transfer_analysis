# The charge_transfer_analysis package

These are scripts to fit the serial CTE signal in the overscan regions
based on a two component model of the trailed charge:

```
q_n = q_0*(Q - b)*exp(-n/tau) + (Q - b)*N*cti**n*(1 - cti) + b
```

Here `n` refers to the n-th overscan column, `q_n` is the mean DN/pixel in
that column, `Q` is the mean DN/pixel in the last imaging column, `b`
is the bias level, and `N` is the number of serial transfers before
the first overscan column (i.e., number of prescan + imaging pixels in
serial direction).  The fit parameters are `q_0`, `tau` and `cti`.

The second term in this expression is the signal expected from the
normal CTI, while the first term is a guess at the behavior of an
additional source of signal smearing which may be due to a "trap" that
releases charge on some time scale `tau`.

This code is intended to be run on high and low flux superflats that
are generated as described in LCA-10103.  Since the bias levels don't
seem to be consistent between those datasets, the `b` value is simply
determined from the final overscan column in each medianed superflat image.

## Set-up

To run this code, one needs the LSST Stack setup and the [eotest
package](https://github.com/lsst-camera-dh/eotest).  At SLAC, where
the Stack is already installed, the download and set up commands from
a rhel6-64 machine running in bash would look like
```
$ scl enable devtoolset-3 bash
$ source /nfs/farm/g/lsst/u1/software/redhat6-x86_64-64bit-gcc44/DMstack/v12_0/loadLSST.bash
$ setup lsst_apps
$ git clone git@github.com:lsst-camera-dh/eotest.git
$ cd eotest
$ setup -r . -j
$ scons opt=3
$ cd .
$ git clone git@github.com:jchiang87/charge_transfer_analysis.git
$ setup -r charge_transfer_analysis -j
$
```
Once the `eotest` and `charge_transfer_analysis` code has been downloaded and built, the
`git clone ...` and `scons ...` steps can be skipped.  The syntax of the setup command
used above is
```
setup -r <directory of package> -j
```
so instead of the relative directory path, as shown in the example, absolute paths can be used.

## People
* [Jim Chiang](https://github.com/jchiang87/charge-transfer-analysis/issues/new?body=@jchiang87) (SLAC)

## License, etc.

This is open source software, available under the BSD license. If you are interested in this project, please do drop us a line via the hyperlinked contact names above, or by [writing us an issue](https://github.com/DarkEnergyScienceCollaboration/charge-transfer-analysis/issues/new).
