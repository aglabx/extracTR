#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 20.02.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com


def get_revcomp(sequence: str) -> str:
    complement = str.maketrans('ATCGNatcgn~[]', 'TAGCNtagcn~][')
    return sequence.translate(complement)[::-1]