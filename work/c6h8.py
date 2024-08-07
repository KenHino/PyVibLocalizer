from pyvib import simulator, read_minfo, write_minfo #geometry, freq, displacement vector


#geom, freq, disp = read_minfo("minfo_sample/hexatriene.minfo", use_trans=True, use_rot=True)
geom, freq, disp = read_minfo("minfo_sample/hexatriene.minfo")

sim_low_freq = simulator(geom, freq[:17], disp[:17])

sim = simulator(geom, freq[17:], disp[17:])

#sim.visualize()
sim.localize(option='Pipek-Mezy', window= 200)
#sim.localize(option='Boys', window=200)
sim.visualize()
write_minfo("hexatriene.out", freq[:17] + sim.freq, disp[:17] + sim.disp)
