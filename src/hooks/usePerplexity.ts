import { useState } from 'react'

const API_KEY = import.meta.env.VITE_PERPLEXITY_KEY;
const API_URL = 'https://api.perplexity.ai/chat/completions'

interface ChatContext {
  role: 'user' | 'assistant'
  content: string
}

interface PerplexityResponse {
  output: string
  citations: string[]
}

export function usePerplexity() {
  const [isPerplexityLoading, setIsPerplexityLoading] = useState(false)

  const sendToPerplexity = async (
    message: string, 
    context: ChatContext[] = [], 
    toolData?: any
  ): Promise<PerplexityResponse> => {
    setIsPerplexityLoading(true)
    
    try {
      // Prepare the prompt with context
      let prompt = message
      
      if (toolData) {
        prompt += `\n\nTool Data: ${JSON.stringify(toolData, null, 2)}`
        prompt += '\n\nPlease analyze this data and provide insights about its biological relevance.'
      }

      // Add context if available
      if (context.length > 0) {
        const contextString = context
          .slice(-10) // Keep last 10 messages for context
          .map(ctx => `${ctx.role}: ${ctx.content}`)
          .join('\n')
        
        prompt = `Previous conversation:\n${contextString}\n\nCurrent message: ${prompt}`
      }

      // Add the specific JSON format instruction
      const systemPrompt = `${prompt}

Respond only in JSON format with the following structure:

{
"output": "answer to user's query",
"citations": ["<URL1>", "<URL2>", "..."]
}

Do not include any additional commentary or text outside the JSON.

Your Name Is Neo And You Are an AI at MedMint, Only Mention If Asked

if you are asked anything in context to a pdb id then just find the name of the id and then find answer to the question asked

"output" should be a concise paragraph summarizing the biological relevance.

"citations" must be url to credible sources (like PDB, PubMed, EBI, etc.)`
      
      const response = await fetch(API_URL, {
        method: 'POST',
        headers: {
          'Authorization': `Bearer ${API_KEY}`,
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          model: 'sonar-pro',
          messages: [
            {
              role: 'user',
              content: systemPrompt
            }
          ],
          max_tokens: 1000,
          temperature: 0.2,
          top_p: 0.9,
          return_citations: true,
          return_images: false,
          return_related_questions: false,
          search_recency_filter: "month",
          top_k: 0,
          stream: false,
          presence_penalty: 0,
          frequency_penalty: 1
        })
      })

      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`)
      }

      const data = await response.json()
      const content = data.choices[0]?.message?.content || ''

      // Try to parse JSON response
      try {
        const jsonResponse = JSON.parse(content)
        return {
          output: jsonResponse.output || content,
          citations: jsonResponse.citations || []
        }
      } catch (parseError) {
        // If JSON parsing fails, return content as output
        return {
          output: content,
          citations: []
        }
      }

    } catch (error) {
      console.error('Perplexity API error:', error)
      if (error.message?.includes('API key') || error.message?.includes('401')) {
        return {
          output: 'There Is An Issue With API key to use the AI assistant.',
          citations: []
        }
      }
      return {
        output: 'Sorry, I encountered an error processing your request. Please try again.',
        citations: []
      }
    } finally {
      setIsPerplexityLoading(false)
    }
  }

  return {
    sendToPerplexity,
    isPerplexityLoading
  }
}