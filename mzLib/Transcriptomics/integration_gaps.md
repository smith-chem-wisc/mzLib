
4. Important integration gap: MetaMorpheus does not consume the generalized property yet.  
   MzlibExtensions.DigestionAgentName() still special-cases only DigestionParams and falls back to the effective agent for RNA. RNA non-specific searches would therefore report singleN/singleC rather than the named RNase until MetaMorpheus is updated.