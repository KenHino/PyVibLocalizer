from pyvib import simulator, read_minfo, write_minfo #geometry, freq, displacement vector


#geom, freq, disp = read_minfo("minfo_sample/wat3.minfo", use_rot=True, use_trans=True)
geom, freq, disp = read_minfo("minfo_sample/wat3.minfo")

sim = simulator(geom, freq, disp)

#sim.visualize(arrow_scale=10)
sim.localize(option='Boys')
#sim.visualize(arrow_scale=10)
#sim.localize(option='Pipek-Mezy')
sim.visualize()

#write_minfo('wat3.out', sim.freq, sim.disp)
