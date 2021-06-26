from sec_groups.secgroups import SecureGroup, repeat_public_base_secret_output


class ElGamal(object):
    """ ElGamal encryption object

    """

    def __init__(self, c1, c2):
        self.c1 = c1
        self.c2 = c2


def encrypt(msgP, k, G, pubkey):
    """Encrypt using ElGamal.

    From [Koblitz1987]: "We further let G in E be a fixed and publicly
    known point. The receiver B chooses a randomly and publishes the key
    aG, while keeping a itself secret. To transmit a message m to B,
    user A chooses a random integer k and sends the pair of points
    (kG, P_m + k(aG))." (Here, P_m is an Encode message m on curve (case: with on EC.))

    """
    return ElGamal(G ^ k, msgP @ (pubkey ^ k))


async def threshold_decrypt(elg_P, sec_x):
    """ Threshold decrypt ElGamal encryption elg_P of P.

    Returns:
    *  SecureGroupElement: Shamir shares of (representation of) group element P.
    """
    # inv_c1_x = elg_P.c1 ^ -sec_x
    group = type(elg_P.c1)
    if hasattr(group, "arithm"):
        sec_grp = SecureGroup(
            group, group.arithm
        )  # Note: For elliptic curves, passing secure arithmetic as second parameter is mandatory.
    else:
        sec_grp = SecureGroup(
            group
        )  # Note: For elliptic curves, passing secure arithmetic as second parameter is mandatory.
    inv_c1_x = repeat_public_base_secret_output(elg_P.c1, -sec_x, sec_grp)
    return inv_c1_x @ elg_P.c2
