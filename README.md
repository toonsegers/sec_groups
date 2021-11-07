[![Build Status](https://app.travis-ci.com/toonsegers/sec_groups.svg)](https://app.travis-ci.com/toonsegers/sec_groups)

# Secure Groups
The 'Secure Groups' Python package implements the Secure Group scheme for several mathematical groups.

## Secure Group scheme
The Secure Group scheme implements finite groups as oblivious data structures. For a given group, the scheme defines the oblivious representation of, and oblivious operations on group elements. Operations include the group law, exponentiation and inversion, random sampling and encoding/decoding.

The oblivious operations are defined by a set of secure multiparty computation (MPC) protocols. Practical protocols are presented for the group of quadratic residues, elliptic curves groups and class groups of imaginary quadratic orders. We demonstrate these protocols in a standard setting for information theoretically secure MPC, tolerating a dishonest minority of passively corrupt parties. 

For our implementation we use the MPyC framework: [https://github.com/lschoe/].

Please find the write-up of this work in Chapter 4 of this Horizon2020 deliverable (published on June 30, 2021): [https://media.voog.com/0000/0042/1115/files/D2.4%20%E2%80%93%20Revision%20of%20Privacy-Enhancing%20Cryptographic%20Primitives%20for%20Ledgers.pdf] 

*Note: This implementation is work-in-progress. (Expect many bugs/issues.) A future version of Secure Groups will be incorporated in MPyC. Until then, please install this work-in-progress version by following the steps explained under Installation.*

### Constant-time extended gcd algorithm and protocol
We introduce a practical protocol to calculate the extended gcd (xgcd) of two secret-shared integers adapting recent work by Bernstein and Yang [https://eprint.iacr.org/2019/266] from the p-adic setting to the finite field setting. This xgcd MPC protocol is a first and of independent interest. We apply it to implement the class group operation in MPC.  

### Conversion protocol for ciphertexts to secret shares
To demonstrate the application of secure groups, we extend a classical threshold cryptosystem with a protocol to convert ciphertexts to secret shares. This functionality enables in- and output to a multiparty computation by communicating one ciphertext over an insecure channel. 

## Installation

This implementation depends on MPyC (version 0.74 or above) and gmpy2.

Install latest version of MPyC:

	git clone https://github.com/lschoe/mpyc
	cd mpyc
	python setup.py install

Install 'gmpy2':

	pip install gmpy2   				# for Linux (first running `apt install libmpc-dev` may be necessary)
	pip install gmpy2-[version etc].whl	# for Windows, see Gohlke's unofficial binaries [https://www.lfd.uci.edu/~gohlke/pythonlibs/]

## Demos

The following demos are included:

* `demo_basic_examples.py` to see examples of different groups (Elliptic curve groups, QR groups, Class groups, etc.);
* `demo_sec_gcd.py` to compute the extended gcd of two (secret shared) integers in constant time;
* `demo_conversion_ed25519.py` to convert ElGamal encryptions to Shamir shares, using the Ed25519 curve group;
* `demo_conversion_qr.py` to convert ElGamal encryptions to Shamir shares, using subgroup of quadratic residues of 2048-bit prime;
* `demo_rubiks.py` to sample random elements from a Rubik's Cube group;

Run the demos as follows:

	cd demos
	python demo_basic_examples.py

## Testing

Run the following commands:

	python -m unittest discover .

## Acknowledgements

This work has received funding from the European Union's Horizon 2020 research and innovation program under grant agreements No 780477 (PRIViLEDGE).
