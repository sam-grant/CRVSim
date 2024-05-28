// TCut goodtrk("trk.status>0");
// TCut CRV1("crvhit[0].sectorType==1");
// TCut KLCRV1("trkfit.sid==200&&trkfit.sindex==1");
// TCut bestfit("trkkl.z0err<1&&trkkl.d0err<1&&trkkl.thetaerr<0.0004&&trkkl.phi0err<0.001&&trk.ndof>=10&&trk.fitcon>0.1&&trk.nactive/trk.nhits>0.99&&trk.nplanes>=4&&trk.nnullambig/trk.nhits<0.2");
// TCut goodCRV("crvhit[0].nHits>0");
// TCut noCRV("crvsummary.nHitCounters==0");
// TCut L1Fiducial("abs(trkfit.pos.X())<2500 && fabs(trkfit.pos.Z()+500)<1500");

TCut goodtrk("kl.status>0");
TCut CRV1("crvcoincs.sectorType==1");
TCut KLCRV1("klfit.sid==200&&klfit.sindex==1");
TCut bestfit("klkl.z0err<1&&klkl.d0err<1&&klkl.thetaerr<0.0004&&klkl.phi0err<0.001&&kl.ndof>=10&&kl.fitcon>0.1&&kl.nactive/kl.nhits>0.99&&kl.nplanes>=4&&kl.nnullambig/kl.nhits<0.2");
TCut goodCRV("crvcoincs.nHits>0");
TCut noCRV("crvsummary.nHitCounters==0");
TCut L1Fiducial("abs(klfit.pos.X())<2500 && fabs(klfit.pos.Z()+500)<1500");