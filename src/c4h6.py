from pyvib import simulator, read_minfo, write_minfo #geometry, freq, displacement vector


#geom, freq, disp = read_minfo("example/butadiene.minfo", use_trans=True, use_rot=True)
geom, freq, disp = read_minfo("example/butadiene.minfo")


sim = simulator(geom, freq, disp)

#sim.visualize()
#sim.localize(option='Pipek-Mezy')
sim.localize(option='Boys')
sim.visualize()
#write_minfo("butadiene.out", sim.freq, sim.disp)
