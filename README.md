## The Story Behind MedMint

### Note- Setup research suite separately, and run it separately!!! use its endpoints in the research directory of this folder, Steps to run the research separately is in the readme file inside the research directory

## Inspiration

**The Problem**: 4,000 diseases have zero FDA-approved treatments. 10 million people die annually from R&D-uncovered diseases. Drug discovery takes 10–15 years, costs $2.6B+, and succeeds only ~12% of the time. Meanwhile, drug monopolies keep life-saving treatments out of reach for patients who need them most.

**The Insight**: AI can generate drug candidates in hours, but most are chemically invalid or fail ADME checks. Worse, when researchers do find something promising, it gets locked in closed silos or owned by large pharma, keeping small labs and independent researchers powerless.

**Our Vision**: What if we could compress years of discovery into weeks, validate every compound rigorously, and register breakthroughs as protected IP on-chain so small pharma and researchers could safely collaborate, keep ownership, and break the monopoly? That's MedMint.
## What it does
MedMint is an AI‑powered breakthrough engine that compresses early drug discovery from 4–5 years to 4–8 weeks.

 - **Compound generation** : Uses TamGen (chemical language model trained on 120M+ compounds) to generate drug candidates for any disease, including rare and neglected ones.

- **Rigorous validation** : Filters candidates with RDKit structural checks, ADME/ADMET profiling, and binding‑affinity prediction so every shortlist is scientifically defensible and patent‑ready.

- **Unified workspace** : A browser‑based platform for molecular design, binding prediction, sequence decoding, structure visualization, chat, whiteboards, shared reports, and task management.

- **On‑chain IP registration** : Every MedMint report is registered as a verifiable IP asset on Story with clear licensing, provenance, and optional royalty tracking.

- **Fair IP layer** : Small pharma and independent researchers can prove priority, license discoveries fairly with automated royalty splits, and partner with big pharma without losing control or credit.

**Impact** : Labs can explore thousands of molecules overnight, get back a rigorous shortlist, and have that work protected and shareable—accelerating and democratizing drug development.

## How we built it
Frontend: React + Next.js for a responsive, real‑time collaborative UI.

Backend: Python‑based services orchestrating AI models and compound pipelines.

AI/ML engine:

TamGen for molecular generation.

RDKit for structural validation and property calculation.

ADME/ADMET models (SwissADME, DeepPurpose) for drug‑likeness profiling.

Binding‑affinity predictors (PDBbind‑trained, Autodock‑Vina) for target validation.

Neo‑style assistant for citation‑backed insights.

Data layer: Supabase for auth, Postgres for reports/metadata, ArangoDB for vector embeddings and similarity search.

Web3/IP layer: Story SDK for registering IP assets, selecting PIL license templates, and wiring in royalty policies.

Deployment: Serverless cloud (e.g., Cloud Run/Lambda) for automatic scaling.

Target or disease → CLM generates candidates → RDKit + ADME filters → binding‑affinity ranking → structured report → one‑click on‑chain IP registration with selected license → optional royalty and partner‑sharing flows.

# MedMint - Decentralized Medical Research Platform

A modern medical research collaboration platform with **Story Protocol** integration for intellectual property (IP) management, licensing, and royalty distribution on the blockchain.

![Story Protocol](https://img.shields.io/badge/Story%20Protocol-Aeneid%20Testnet-purple)
![React](https://img.shields.io/badge/React-18.3-blue)
![TypeScript](https://img.shields.io/badge/TypeScript-5.5-blue)
![Vite](https://img.shields.io/badge/Vite-5.4-yellow)

## Quick Start

### Prerequisites

- **Node.js** 18+ 
- **npm** or **yarn**
- **MetaMask** or compatible Web3 wallet
- **Story Protocol Aeneid testnet** configured in wallet

### Installation

```bash
# Clone the repository
git clone <your-repo-url>
cd baitmint

# Install dependencies
npm install

# Create environment file
cp .env.example .env
```

### Environment Variables

Create a `.env` file in the root directory:

```env
# Supabase Configuration
VITE_SUPABASE_URL=your_supabase_url
VITE_SUPABASE_ANON_KEY=your_supabase_anon_key

# Pinata IPFS (for metadata storage)
VITE_PINATA_JWT=your_pinata_jwt_token

# WalletConnect Project ID (optional)
VITE_WALLETCONNECT_PROJECT_ID=your_walletconnect_project_id
```

### Running the App

```bash
# Development mode
npm run dev

# Build for production
npm run build

# Preview production build
npm run preview
```

The app will be available at `http://localhost:5173`

---

## Story Protocol Integration

MedMint integrates with **Story Protocol** on the **Aeneid Testnet** to provide decentralized IP management for medical research documents.

### Network Configuration

| Property | Value |
|----------|-------|
| **Network** | Story Aeneid Testnet |
| **Chain ID** | 1315 |
| **RPC URL** | https://aeneid.storyrpc.io |
| **Explorer** | https://aeneid.storyscan.xyz |
| **Currency** | $IP |

### Contract Addresses (Aeneid Testnet)

| Contract | Address |
|----------|---------|
| **SPG NFT Contract** | `0xc32A8a0FF3beDDDa58393d022aF433e78739FAbc` |
| **WIP Token** | `0x1514000000000000000000000000000000000000` |
| **Royalty Policy LAP** | `0xBe54FB168b3c982b7AaE60dB6CF75Bd8447b390E` |

---

## Story Protocol SDK Features Used

### 1. **StoryClient Initialization**
```typescript
import { StoryClient, StoryConfig } from '@story-protocol/core-sdk'

const config: StoryConfig = {
  account: walletClient.account,
  transport: custom(walletClient.transport),
  chainId: 'aeneid',
}
const client = StoryClient.newClient(config)
```

### 2. **IP Asset Registration** (`client.ipAsset.registerIpAsset`)

Registers research reports as IP Assets on Story Protocol with metadata stored on IPFS.

**Features:**
- Mints an SPG NFT representing the IP
- Uploads metadata to IPFS via Pinata
- Attaches license terms during registration
- Computes SHA-256 hash of metadata for verification

```typescript
const response = await client.ipAsset.registerIpAsset({
  nft: {
    type: 'mint',
    spgNftContract: SPG_NFT_CONTRACT,
  },
  ipMetadata: {
    ipMetadataURI,
    ipMetadataHash,
    nftMetadataURI,
    nftMetadataHash,
  },
  licenseTermsData: [...], // Optional license terms
})
```

also in .env file add these 
```
VITE_SUPABASE_URL=https://
VITE_SUPABASE_ANON_KEY= Anon Key Here
VITE_PERPLEXITY_KEY = Sonar API Key Here
VITE_PINATA_API_SECRET=
VITE_PINATA_API_KEY=
VITE_PINATA_JWT=
```
### 3. **License Terms Configuration**

MedMint supports three license types with Programmable IP License (PIL) terms:

| License Type | Commercial Use | Derivatives | Minting Fee | Rev Share |
|--------------|---------------|-------------|-------------|-----------|
| **Open Research** | ❌ No | ✅ Yes | Free | 0% |
| **Commercial** | ✅ Yes | ❌ No | 1 $IP | 0% |
| **Commercial Remix** | ✅ Yes | ✅ Yes | 0.5 $IP | 5% |

**PIL Terms Structure:**
```typescript
{
  transferable: true,
  royaltyPolicy: ROYALTY_POLICY_LAP,
  defaultMintingFee: parseEther('1'), // 1 $IP
  commercialUse: true,
  commercialAttribution: true,
  derivativesAllowed: false,
  currency: WIP_TOKEN_ADDRESS,
  // ... more terms
}
```

### 4. **License Token Minting** (`client.license.mintLicenseTokens`)

Allows others to mint license tokens to use the IP. Minting fee goes to the IP's royalty vault.

```typescript
const response = await client.license.mintLicenseTokens({
  licenseTermsId: '1',
  licensorIpId: ipId,
  receiver: receiverAddress,
  amount: 1,
})
```

### 5. **Royalty Management**

#### Get Royalty Vault Address (`client.royalty.getRoyaltyVaultAddress`)
```typescript
const vaultAddress = await client.royalty.getRoyaltyVaultAddress(ipId)
```

#### Check Claimable Revenue (`client.royalty.claimableRevenue`)
```typescript
const claimable = await client.royalty.claimableRevenue({
  ipId: ipId,
  claimer: address,
  token: WIP_TOKEN_ADDRESS,
})
```

#### Claim All Revenue (`client.royalty.claimAllRevenue`)
```typescript
const response = await client.royalty.claimAllRevenue({
  ancestorIpId: ipId,
  claimer: address,
  childIpIds: [],
  royaltyPolicies: [ROYALTY_POLICY_LAP],
  currencyTokens: [WIP_TOKEN_ADDRESS],
})
```

### 6. **Vault Snapshot** (Direct Contract Call)

For funds that haven't been snapshotted yet, we call the vault's `snapshot()` function directly:

```typescript
const snapshotTxHash = await walletClient.writeContract({
  address: vaultAddress,
  abi: [{ name: 'snapshot', type: 'function', ... }],
  functionName: 'snapshot',
})
```

---

## Architecture

### Core Components

```
src/
├── hooks/
│   └── useStoryProtocol.ts    # Story Protocol integration hook
├── components/
│   ├── ReportsPage.tsx        # Reports management with IP features
│   ├── RoyaltyPanel.tsx       # Royalty viewing and claiming UI
│   └── Web3Provider.tsx       # Wagmi/Web3 provider setup
└── lib/
    └── web3.ts                # Web3 configuration
```

### Hook: `useStoryProtocol`

The main hook exports the following functions:

| Function | Description |
|----------|-------------|
| `registerIP(metadata)` | Register a report as IP with optional license terms |
| `mintLicenseToken(ipId, termsId, receiver, amount)` | Mint license tokens for an IP |
| `getRoyaltyInfo(ipId)` | Get vault balance, claimable amount, and royalty details |
| `claimRoyalties(ipId)` | Claim available royalties |
| `snapshotAndClaim(ipId)` | Snapshot vault and claim funds |
| `getVaultBalance(ipId)` | Debug function to check vault vs claimable |

### Exported Types

```typescript
export type LicenseType = 'none' | 'open' | 'commercial' | 'commercial-remix'

export interface RoyaltyInfo {
  ipId: string
  claimableAmount: string
  vaultBalance: string
  unsnapshotedAmount: string
  hasUnsnapshotedFunds: boolean
  royaltyVaultAddress?: string
  // ...
}

export interface RegisterIPResult {
  txHash: string
  ipId: string
  tokenId: string
  licenseTermsId?: string
}
```

---

## Permission-Based Access

MedMint implements role-based access control for Story Protocol features:

| Action | Lab Admin/Owner | Lab Member |
|--------|-----------------|------------|
| Connect Wallet | ✅ | ✅ |
| View Reports | ✅ | ✅ |
| Register as IP | ✅ | ❌ |
| Share IP (Mint License) | ✅ | ✅ |
| View Royalties | ✅ | ❌ |
| Claim Royalties | ✅ | ❌ |

---

## User Flows

### 1. Register Research as IP

1. Navigate to **Reports** page
2. Click on a report to view details
3. Click **"Register as IP"** (admin only)
4. Select license type:
   - No License
   - Open Research (Non-Commercial)
   - Commercial Use (1 $IP minting fee)
   - Commercial Remix (0.5 $IP + 5% rev share)
5. Confirm transaction in wallet
6. View IP on Story Explorer

### 2. Share IP / Mint License

1. Click **"Share IP (Mint License)"** on a registered report
2. Enter recipient wallet address
3. Confirm transaction (recipient pays minting fee)
4. License token is minted to recipient

### 3. Claim Royalties

1. Click **"View Royalties"** on a registered report (admin only)
2. View vault balance and claimable amount
3. If funds are unsnapshotted, click **"Snapshot & Claim"**
4. Otherwise, click **"Claim"**
5. Royalties are transferred to your wallet

---

## Tech Stack

- **Frontend**: React 18, TypeScript, Vite
- **Styling**: Tailwind CSS
- **Web3**: wagmi, viem
- **Story Protocol**: @story-protocol/core-sdk v1.4.2
- **Storage**: IPFS (via Pinata), Supabase
- **Icons**: Lucide React

---

## License

MIT License

---

## Resources

- [Story Protocol Documentation](https://docs.story.foundation/)
- [Story Protocol SDK](https://github.com/storyprotocol/sdk)
- [Aeneid Testnet Explorer](https://aeneid.explorer.story.foundation/)
- [Aeneid Faucet](https://faucet.story.foundation/)


## Challenges we ran into

Speed vs accuracy: Generating thousands of molecules is easy; returning a high‑quality shortlist fast is hard. Multi‑stage filtering, GPU batching, and caching kept latency low without sacrificing rigor.

Complex IP plumbing, simple UX: Story’s licensing and royalty modules are powerful but complex. MedMint hides that behind clear license presets and human‑readable summaries.

Noisy biological data: ADME datasets are incomplete and inconsistent. Cross‑model validation and similarity‑based sanity checks improved reliability.

Real‑time collaboration: Keeping chats, whiteboards, and reports in sync under load required careful state management and conflict‑free data structures.

Provenance mapping: Ensuring every report version links correctly to its on‑chain IP asset meant designing explicit versioning and transaction tracking.

## Accomplishments that we're proud of
100x faster early discovery: From target to vetted shortlist in minutes to hours instead of months to years.

Scientifically grounded results: Shortlists validated against peer‑reviewed ADME and binding‑affinity models, not just black‑box scores.

IP that actually works for the little guy: Small labs can now register discoveries, set licenses, and negotiate fairly with big pharma.

One workspace, many tools replaced: MedMint subsumes docking tools, drawing tools, spreadsheets, and chat into a single coordinated environment.

Real‑world validation: Early feedback and interest from pharma researchers and analysts who stress‑tested the outputs and found them credible.

## What we learned

Owning IP changes behavior: When small teams know their work is provably theirs, they are more willing to share, collaborate, and publish intermediate results.

Trust beats pure speed: Scientists accept a slightly slower model if they understand the validation pipeline and can cite underlying literature.

Rare diseases need better tools: The biggest excitement came from teams working on conditions big pharma often ignores.

UX for Web3 must feel Web2: Hiding wallets, gas, and protocol jargon behind clean flows is essential for adoption in life sciences.

Collaboration is a feature, not a nice‑to‑have: Shared context (chat + whiteboards + reports) is as important as the model itself.

## What's next for MedMint

Royalty and revenue dashboards: Show labs what they can claim and what flows from licensed or derivative use of their IP.

Partner licensing flows: One‑click “share this discovery with partner X” via license tokens and templated agreements.

More powerful AI agents: Autonomous workflows that iteratively generate, filter, and optimize compounds for specific constraints.

Wet‑lab integration: Connect in‑silico runs with lab automation and sensor data for closed‑loop optimization.

Discovery DAO and marketplace: Pool results from many labs under shared governance and let anyone license data, models, or validated hits—reducing monopolies and making cures more accessible worldwide.


Ethics and Safety: We plan to consult with pharmacologists, ethicists, and regulatory advisors to ensure the platform supports responsible discovery and transparent AI decision-making.
