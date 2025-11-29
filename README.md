## The Story Behind MedMint

### Note- Setup research suite separately, and run it separately!!! use its endpoints in the research directory of this folder

It started with a question we couldn’t shake off:
“What if the drug that could save someone’s life hasn’t been discovered yet, not because it’s impossible, but because the right compound was never guessed correctly by hit and trial, the right target was never predicted, or the compound never passed clinical trials?”

We spoke to my uncle, who is chief radiologist at Apollo, one of the best hospital chains in India, from students in academia & to scientists in early-stage biotech — and all echoed the same struggle:
The science was there. The ideas were bold. But the process? Painfully slow and high error rates.

**Behind every stalled experiment is a family, a person, someone's parents, someone's children waiting hopefully a cure would be there for them. Slow research doesn’t just delay science — it delays hope. Families uproot their lives chasing second opinions, while potential breakthroughs sit trapped in outdated systems and scattered workflows.**

In India, families waiting for enzyme replacement therapy for rare genetic disorders have faced months-long delays, even after funding was approved and in at least 10 tragic cases, children died before treatment could begin [source](https://www.theippress.com/2024/10/13/justice-for-the-rare-unraveling-delays-access-and-the-fight-for-life-saving-drugs-in-india/)

#How MedMint Solve This Major Issue
MedMint stands for - Medicinal & Molecular Intelligence Toolkit

In a sea of quick answers and chatbots, MedMint chose the hard path: tackling one of the most complex, expensive, and uncertain stages in healthcare: early-stage drug discovery.

We didn’t just build a demo. We built an AI-powered pipeline that takes you from a disease target to real medicinal compound generation and explainable evaluation, all in one platform. Using models like TamGen, chemBERTa, and DeepPurpose, we made scientific discovery searchable, explainable, and actionable.

**This isn’t just about tech — it’s about saving years of research, millions in cost, and potentially giving hope to patients who currently have no cure.**

MedMint is a seamless research and collaboration suite designed to accelerate early-stage drug discovery.

It helps researchers generate medicinal compounds, predict targets, analyze binding affinities, decode protein sequences, and visualize molecular structures in 3D, all in one unified workspace,  offering real-time chat, collaborative whiteboards, shared task systems, and a lab-wide report repository. It ensures that researchers, no matter where they are, can work as if they’re side by side.

MedMint streamlines the most resource-heavy and uncertain stages of early drug discovery — from identifying a viable target to generating and evaluating candidate compounds. This drastically reduces time, cost, and cognitive overhead for researchers, bridging gaps between silos and boosting the chances of clinical success.

## What it does
Drug discovery takes around 10-15 years [Source: DiMasi, J. A., Grabowski, H. G., & Hansen, R. W. (2016). Innovation in the pharmaceutical industry: New estimates of R&D costs. Journal of Health Economics, 47, 20–33.] MedMint accelerates drug discovery by unifying the entire pipeline — from molecule creation to deep biological insight — into a single, intuitive system. It’s not just faster — it’s smarter, more connected, and built for the urgency of real-world healthcare, saving 3-5 years of development time and plans to save more in future iterations

At the core of medicinal compound generation is TamGen, a cutting-edge generative AI model trained on over 120 million bioactive molecules. Unlike traditional methods, TamGen generates target-specific, synthesizable compounds — and has already produced novel molecules that demonstrated strong biological activity in real-world validations [Source: nature article](https://www.nature.com/articles/s41467-024-53632-4). In a research on TamGe, TamGen predicted and identified novel compounds that act as inhibitors of the ClpP protease from Mycobacterium tuberculosis, the bacterium responsible for tuberculosis (TB), the identified compounds include novel scaffolds such as benzenesulfonamide and diphenylurea groups, which differ significantly from existing TB protease inhibitors like Bortezomib, potentially offering improved bioavailability and molecular stability.

By combining this with built-in ADME profiling, binding affinity prediction, and explainable AI tools, MedMint helps researchers prioritize compounds with true therapeutic potential — cutting 3–5 years from the timeline and millions of dollars without compromising scientific rigor & and hopefully making medicines cheaper

[Link To MedMint Presentation](https://www.canva.com/design/DAGrw_xpBVk/W7Z_QLLCTwrU5w66jy1HsQ/view?utm_content=DAGrw_xpBVk&utm_campaign=designshare&utm_medium=link2&utm_source=uniquelinks&utlId=h1a34c55d90)

## Tech We Used
- Supabase (For Database, Edge Functions, Triggers, Storage Bucket, RLS Policies, Broadcast, Auth, Email Confirmation Of Auth)

- Entri (For Domain & Domain Management)

- Netlify for deployment

At the heart of MedMint it is built on Bolt, which served as the foundation for our entire platform — powering real-time collaboration, user management, task systems, and secure lab environments. Every core interaction — from report sharing to whiteboard syncing and team messaging — runs on Bolt’s robust infrastructure, enabling fast, secure, and scalable performance right out of the box. To extend Bolt’s capabilities into advanced drug discovery, we integrated specialized AI-powered microservice, These microservices were modularly deployed on google cloud and seamlessly connected via Bolt, making our system not just powerful, but incredibly flexible. Bolt helped us create an event-driven architecture allowed us to orchestrate tasks between collaborative lab tools and scientific engines — ensuring everything from a generated molecule to its ADME profile was tied back into the research flow in real time. Bolt’s flexibility let us route data, trigger events, and deliver responses inside the collaborative lab flow, without breaking a sweat.

If MedMint is the lab of the future, then **Bolt is the electricity that powers every single experiment**.

## Challenges we ran into
Building MedMint in under a month meant pushing every tool to its limits — and sometimes beyond

Bolt.new’s project size limits: As our application scaled and the features grew, we hit Bolt’s project size limits. We had to optimize data handling and offload some large state interactions to external stores.

Prompt engineering: While integrating AI for research queries and compound generation, we realized how critical prompt design is. It took iterative tuning to ensure domain-specific queries returned accurate, trustworthy, and citable responses.

Desync between Bolt and server state: During late-night sprints, we built some of our most innovative features — only to find they had disappeared the next day. Turns out, syncing issues between Bolt’s editor and server state led to features vanishing silently, without persistent saves. We had to fork from our last stable backup and reimplement them from scratch, losing precious hours but gaining hard-earned lessons on state durability and backup hygiene.

Supabase Realtime: We were using supabase realtime for real time chat messages but we wondered why it was not working, turns out that feature is not available publically to everyone yet

## Accomplishments that we're proud of
MedMint is a crown full of jewels, the platform itself is a very proud thing to have built, but the things we are the most proud of are-

Designed a platform with real-world impact: something that could save researchers years and bring patients hope sooner — and got early validation from medical professionals and domain experts.

Enabled real-time collaboration — shared whiteboards, chat, task tracking, and a shared report repository — all functioning like a digital research lab.

Integrated advanced drug discovery tools (compound generation, target prediction, ADME profiling, amino acid decoding, binding affinity) as modular microservices — and made them feel like one connected platform.

Built a full-stack research suite from scratch in under a month— with AI, real-time collaboration, molecular visualization, and report generation all working together seamlessly.

Used TamGen to generate target-aware compounds with drug-like properties, backed by real scientific research and explainable outputs.

And new users signing up already

## What we learned

AI can’t replace researchers — but it can empower them. We learned that true innovation comes when AI is used to assist, not replace. MedMint isn't about automating science, but about giving researchers superpowers: faster insights, stronger predictions, and less time lost on the wrong leads.

AI alone isn’t enough. Even powerful models like TamGen need careful inputs, validation, and scientific interpretation. We gained a deep respect for the balance between computation and context.

Tools can break — backups shouldn’t. Our Bolt desync taught us to expect the unexpected. Having fallback systems and version control saved our project (and our sanity).

Drug discovery is more than code. It’s about understanding human urgency, regulatory reality, and the daily friction researchers face. Talking to real scientists gave us clarity on what actually matters.

The future of science is collaborative. What we’re building is not just a research suite — it’s a space where people can come together and make discoveries that no one could make alone.

## What's next for MedMint

A future where researchers can join public or private labs, access datasets, contribute models, and share reports — building a living ecosystem of scientific progress.

Tighter integration with tools like ElevenLabs (for voice), AlphaFold (for structure), and public compound databases to push MedMint from a discovery engine to a decision-making ally in preclinical development.

Ethics and Safety: We plan to consult with pharmacologists, ethicists, and regulatory advisors to ensure the platform supports responsible discovery and transparent AI decision-making.
