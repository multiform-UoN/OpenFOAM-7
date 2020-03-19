#!/bin/bash
cat $1 | grep areaIntegrate | cut -d' ' -f9 | tr -d ',(' > fluxinlet
